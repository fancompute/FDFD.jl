export Modulator
export solve_modulation_TM
export assign_mod_delta!, assign_mod_delta_func!
export assign_mod_phi!, assign_mod_phi_func!

using Pardiso, TimerOutputs

mutable struct Modulator
    geom::Geometry2D
    Omega::Real
    Nsb::Integer
    epsr_delta::Array{Complex,2}
    epsr_delta_phi::Array{Real,2}

    function Modulator(geom::Geometry2D, Omega::Real, Nsb::Integer)
        epsr_delta = zeros(Complex128, geom.N)
        epsr_delta_phi = zeros(Float64, geom.N)
        return new(geom, Omega, Nsb, epsr_delta, epsr_delta_phi)
    end
end

function assign_mod_delta!(mod::Modulator, region, value)
    mask = [region(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    mod.epsr_delta[mask] = value;
end

function assign_mod_delta_func!(mod::Modulator, region, value_func)
    mask = [region(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    value_computed = [value_func(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    mod.epsr_delta[mask] = value_computed[mask];
end

function assign_mod_phi!(mod::Modulator, region, value)
    mask = [region(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    mod.epsr_delta_phi[mask] = value;
end

function assign_mod_phi_func!(mod::Modulator, region, value_func)
    mask = [region(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    value_computed = [value_func(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    mod.epsr_delta_phi[mask] = value_computed[mask];
end

function solve_modulation_TM(mod::Modulator, omega0; verbose=true)
    const to = TimerOutput()

    @timeit to "setup" begin
    println("# Solver: setup...");

    geom = mod.geom;
    M = prod(geom.N);
    Nsb = mod.Nsb

    println("# Solver: ", @sprintf("%.1E",(2*Nsb+1)*M), " unknowns");

    n_sb = -Nsb:1:Nsb; 
    omega = omega0 + mod.Omega*n_sb; 

    Ez = zeros(Complex128, 2*Nsb+1, geom.N[1], geom.N[2]);
    Hx = zeros(Complex128, 2*Nsb+1, geom.N[1], geom.N[2]);
    Hy = zeros(Complex128, 2*Nsb+1, geom.N[1], geom.N[2]);

    T_eps = spdiagm(epsilon0*geom.epsr[:]); #TODO: check reshape vs [:]
    T_delta = spdiagm(epsilon0*mod.epsr_delta[:]); #TODO: check reshape vs [:]
    T_phi = spdiagm(exp.(1im*mod.epsr_delta_phi[:])); #TODO: check reshape vs [:]

    # Construct derivates
    Dxb = dws("x", "b", geom);
    Dxf = dws("x", "f", geom);
    Dyb = dws("y", "b", geom);
    Dyf = dws("y", "f", geom); 

    # Reshape Mz into a vector
    b0 = 1im*omega[Nsb+1]*geom.src[:]; #TODO: check reshape vs [:]
    b = zeros(Complex128, M*(2*Nsb+1),1); 
    b[(Nsb*M)+1:(Nsb+1)*M,1] = b0; 
    end

    println("# Solver: starting matrix assembly");
    @timeit to "matrix assembly" begin
    As  = Array{SparseMatrixCSC}(2*Nsb+1);
    Sxf = Array{SparseMatrixCSC}(2*Nsb+1);
    Sxb = Array{SparseMatrixCSC}(2*Nsb+1);
    Syf = Array{SparseMatrixCSC}(2*Nsb+1);
    Syb = Array{SparseMatrixCSC}(2*Nsb+1);

    println("#         system matrices...");
    for i = 1:(2*Nsb + 1)
        @timeit to "S_i" (Sxf[i], Sxb[i], Syf[i], Syb[i]) = S_create(omega[i], geom.N, geom.Npml, geom.xrange, geom.yrange);
        @timeit to "A_i" As[i] = Sxb[i]*Dxb*mu0^-1*Sxf[i]*Dxf + Syb[i]*Dyb*mu0^-1*Syf[i]*Dyf + omega[i]^2*T_eps;
    end
    if Nsb > 0
        println("#         coupling matrices...");
        @timeit to "C_p" begin
            template_Cp = spdiagm([0.5*omega[1:end-1].^2], [1], 2*Nsb+1, 2*Nsb+1);
            Cp = kron(template_Cp, T_delta*conj(T_phi));
        end
        @timeit to "C_m" begin
            template_Cm = spdiagm([0.5*omega[2:end].^2],  [-1], 2*Nsb+1, 2*Nsb+1);
            Cm = kron(template_Cm, T_delta*T_phi);
        end
        println("#         assembling...");
        @timeit to "A" A = blkdiag(As...) + Cp + Cm;
    else
        A = As[1]; 
    end
    end


    println("# Solver: solving...");
    @timeit to "solve" begin
    pardiso_success = false;
    try
        ez = zeros(Complex128, M*(2*Nsb+1),1); 
        ps = PardisoSolver();
        if verbose
        	set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
        end
        set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM);
        set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON);
        set_solver!(ps, Pardiso.DIRECT_SOLVER);
        pardisoinit(ps);
        solve!(ps, ez, A, b);
        pardiso_success = true;
        println("# Solver: pardiso performed %d iterative refinement steps", get_iparm(ps, 7));
        #println("# Solver: maximum residual for the solution X is %0.3g", maximum(abs(A*ez-b)));
    catch
        println("# Solver: pardiso failed, falling back to lufact()");
    end

    if ~pardiso_success
        ez = lufact(A)\b;
    end
    end

    println("# Solver: post processing...");
    @timeit to "post processing" begin
    for i = 1:(2*Nsb+1)
        Ez[i, :, :] = reshape(ez[(i-1)*M+1:i*M], geom.N);
        Hx[i, :, :] = reshape(-1/1im/omega[i]*mu0^-1*Syf[i]*Dyf*ez[(i-1)*M+1:i*M], geom.N); 
        Hy[i, :, :] = reshape(1/1im/omega[i]*mu0^-1*Sxf[i]*Dxf*ez[(i-1)*M+1:i*M], geom.N); 
    end
    end

    if verbose
    	show(to);
    	println("");
    end
    
    return (Ez, Hx, Hy, omega)
end
