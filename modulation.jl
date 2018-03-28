export Modulator
export solve_modulation_TM
export assign_mod_delta!

using Pardiso, TimerOutputs

mutable struct Modulator
    geom::Geometry2D
    Ω::Real
    nsidebands::Integer
    Δϵᵣ::Array{Complex,2}

    function Modulator(geom::Geometry2D, Ω::Real, nsidebands::Integer)
        Δϵᵣ = zeros(Complex128, geom.N)
        return new(geom, Ω, nsidebands, Δϵᵣ)
    end
end

function assign_mod_delta!(mod::Modulator, region, value)
    mask = [region(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    if iscallable(value)
        value_assigned = [value(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    else
        value_assigned = value
    end
    mod.Δϵᵣ[mask] = value;
end

function solve_modulation_TM(mod::Modulator, ω₀; verbose=true)
    geom = mod.geom;
    nsidebands = mod.nsidebands;
    n = -nsidebands:1:nsidebands; 
    N = geom.N;
    M = prod(N);
    Ω = mod.Ω;
    ω = ω₀ + Ω*n; 

    const to = TimerOutput()
    @timeit to "setup" begin
    println("# Solver: setup...");

    println("# Solver: ", @sprintf("%.2E", (2*nsidebands+1)*M), " unknowns");

    Ez = zeros(Complex128, 2*nsidebands+1, N[1], N[2]);
    Hx = zeros(Complex128, 2*nsidebands+1, N[1], N[2]);
    Hy = zeros(Complex128, 2*nsidebands+1, N[1], N[2]);

    Tϵ = spdiagm(ϵ₀*geom.ϵᵣ[:]); #TODO: check reshape vs [:]
    TΔϵ = spdiagm(ϵ₀*mod.Δϵᵣ[:]); #TODO: check reshape vs [:]

    # Construct derivates
    δxb = δ("x", "b", geom);
    δxf = δ("x", "f", geom);
    δyb = δ("y", "b", geom);
    δyf = δ("y", "f", geom); 

    # Reshape Mz into a vector
    b0 = 1im*ω₀*geom.src[:]; #TODO: check reshape vs [:]
    b = zeros(Complex128, M*(2*nsidebands+1),1); 
    b[(nsidebands*M)+1:(nsidebands+1)*M,1] = b0; 
    end

    println("# Solver: starting matrix assembly");
    @timeit to "matrix assembly" begin
    As  = Array{SparseMatrixCSC}(2*nsidebands+1);
    Sxf = Array{SparseMatrixCSC}(2*nsidebands+1);
    Sxb = Array{SparseMatrixCSC}(2*nsidebands+1);
    Syf = Array{SparseMatrixCSC}(2*nsidebands+1);
    Syb = Array{SparseMatrixCSC}(2*nsidebands+1);

    println("#         system matrices...");
    for i = 1:(2*nsidebands + 1)
        @timeit to "S_i" (Sxf[i], Sxb[i], Syf[i], Syb[i]) = S_create(ω[i], N, geom.Npml, geom.xrange, geom.yrange);
        @timeit to "A_i" As[i] = Sxb[i]*δxb*μ₀^-1*Sxf[i]*δxf + Syb[i]*δyb*μ₀^-1*Syf[i]*δyf + ω[i]^2*Tϵ;
    end
    if nsidebands > 0
        println("#         coupling matrices...");
        @timeit to "C_p" begin
            template_Cp = spdiagm([0.5*ω[1:end-1].^2], [1], 2*nsidebands+1, 2*nsidebands+1);
            Cp = kron(template_Cp, conj(TΔϵ));
        end
        @timeit to "C_m" begin
            template_Cm = spdiagm([0.5*ω[2:end].^2],  [-1], 2*nsidebands+1, 2*nsidebands+1);
            Cm = kron(template_Cm, TΔϵ);
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
        ez = zeros(Complex128, M*(2*nsidebands+1),1); 
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
    for i = 1:(2*nsidebands+1)
        Ez[i, :, :] = reshape(ez[(i-1)*M+1:i*M], N);
        Hx[i, :, :] = reshape(-1/1im/ω[i]*μ₀^-1*Syf[i]*δyf*ez[(i-1)*M+1:i*M], N); 
        Hy[i, :, :] = reshape(1/1im/ω[i]*μ₀^-1*Sxf[i]*δxf*ez[(i-1)*M+1:i*M], N); 
    end
    end

    if verbose
    	show(to);
    	println("");
    end
    
    return (Ez, Hx, Hy, ω)
end
