export solve_modulation_TM
export assign_mod_delta!, assign_mod_delta_func!
export assign_mod_phi!, assign_mod_phi_func!

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

function solve_modulation_TM(mod::Modulator, omega0)
    geom = mod.geom;
    M = prod(geom.N);
    Nsb = mod.Nsb

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

    As  = Array{SparseMatrixCSC}(2*Nsb+1);
    Sxf = Array{SparseMatrixCSC}(2*Nsb+1);
    Sxb = Array{SparseMatrixCSC}(2*Nsb+1);
    Syf = Array{SparseMatrixCSC}(2*Nsb+1);
    Syb = Array{SparseMatrixCSC}(2*Nsb+1);
    for i = 1:(2*Nsb + 1)
        (Sxf[i], Sxb[i], Syf[i], Syb[i]) = S_create(omega[i], geom.N, geom.Npml, geom.xrange, geom.yrange);
        As[i] = Sxb[i]*Dxb*mu0^-1*Sxf[i]*Dxf + Syb[i]*Dyb*mu0^-1*Syf[i]*Dyf + omega[i]^2*T_eps;
    end

    if Nsb > 0
        template_Cp = spdiagm([0.5*omega[1:end-1].^2], [1], 2*Nsb+1, 2*Nsb+1)
        template_Cm = spdiagm([0.5*omega[2:end].^2],  [-1], 2*Nsb+1, 2*Nsb+1)
        Cp = kron(template_Cp, T_delta*conj(T_phi));
        Cm = kron(template_Cm, T_delta*T_phi);
        A = blkdiag(As...)+Cp+Cm;
    else
        A = As[1]; 
    end

    if solver_pardiso
        @time ez = solve(handle_ps, A, b);
    else
        @time ez = lufact(A)\b;
    end

    for i = 1:(2*Nsb+1)
        Ez[i,:,:] = reshape(ez[(i-1)*M+1:i*M], geom.N);
        Hx[i,:,:] = reshape(-1/1im/omega[i]*mu0^-1*Syf[i]*Dyf*ez[(i-1)*M+1:i*M], geom.N); 
        Hy[i,:,:] = reshape(1/1im/omega[i]*mu0^-1*Sxf[i]*Dxf*ez[(i-1)*M+1:i*M], geom.N); 
    end

    return (Ez, Hx, Hy, omega)
end