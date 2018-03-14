export solve_modulation_TM

function solve_modulation_TM(omega0, Omega, Nsb, xrange, yrange, eps_r, mod_reg, mod_phi, Jz0, Npml)
    N = size(eps_r);
    M = prod(N); 

    n_sb = -Nsb:1:Nsb; 

    Ez = zeros(Complex128, 2*Nsb+1, N[1], N[2]);
    Hx = zeros(Complex128, 2*Nsb+1, N[1], N[2]);
    Hy = zeros(Complex128, 2*Nsb+1, N[1], N[2]);
    
    omega = omega0 + Omega*n_sb; 

    T_eps = spdiagm(epsilon0*eps_r[:]); #TODO: check reshape vs [:]
    T_delta = spdiagm(epsilon0*mod_reg[:]); #TODO: check reshape vs [:]
    T_phi = spdiagm(exp.(1im*mod_phi[:])); #TODO: check reshape vs [:]

    # Construct derivates
    Dxb = dws("x", "b", N, xrange, yrange);
    Dxf = dws("x", "f", N, xrange, yrange);
    Dyb = dws("y", "b", N, xrange, yrange);
    Dyf = dws("y", "f", N, xrange, yrange); 

    # Reshape Mz into a vector
    b0 = 1im*omega[Nsb+1]*Jz0[:]; #TODO: check reshape vs [:]
    b = zeros(Complex, M*(2*Nsb+1),1); 
    b[(Nsb*M)+1:(Nsb+1)*M,1] = b0; 

    As = Array{SparseMatrixCSC}(2*Nsb+1);
    for i = 1:(2*Nsb + 1)
        (Sxf, Sxb, Syf, Syb) = S_create(omega[i], N, Npml, xrange, yrange);
        As[i] = Sxb*Dxb*mu0^-1*Sxf*Dxf + Syb*Dyb*mu0^-1*Syf*Dyf + omega[i]^2*T_eps;
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
    ez = A\b;

    for i = 1:(2*Nsb+1)
        Ez[i,:,:] = reshape(ez[(i-1)*M+1:i*M, 1], N);
        (Sxf, Sxb, Syf, Syb) = S_create(omega[i], N, Npml, xrange, yrange);
        Hx[i,:,:] = -1/1im/omega[i]*mu0^-1*Syf*Dyf*ez[(i-1)*M+1:i*M, 1]; 
        Hy[i,:,:] = 1/1im/omega[i]*mu0^-1*Sxf*Dxf*ez[(i-1)*M+1:i*M, 1]; 
    end

    return (Ez, Hx, Hy, omega)
end