function solve_eigen_1D_omega(pol, beta, omega_est, Neigs, xrange, eps_r)
    N = size(eps_r);
    M = prod(N);

    eps_x = grid_average(eps_r, "x");
    T_eps = spdiagm(epsilon0*eps_r[:]);
    T_eps_x = spdiagm(epsilon0*eps_x[:]);

    Dxb = dws("x", "b", (N[1],1), xrange, (-Inf, Inf));
    Dxf = dws("x", "f", (N[1],1), xrange, (-Inf, Inf));

    #if pol == "TM"
    A = mu0^-1*T_eps.^-1*beta^2 + mu0^-1*T_eps.^-1*Dxf*Dxb;
    #elseif pol == "TE"
    #    A = omega^2*mu0*T_eps + T_eps*Dxf*T_eps_x.^-1*Dxb;
    #else
    #    error("Invalid polarization specified!");
    #end

    (omega_sqr, output_vectors) = eigs(A, nev=Neigs, sigma=omega_est^2);
    omega = sqrt.(omega_sqr);
    return (omega, output_vectors)
end

function solve_eigen_1D(pol, omega, beta_est, Neigs, xrange, eps_r)
    N = size(eps_r);
    M = prod(N);

    eps_x = grid_average(eps_r, "x");
    T_eps = spdiagm(epsilon0*eps_r[:]);
    T_eps_x = spdiagm(epsilon0*eps_x[:]);

    Dxb = dws("x", "b", (N[1],1), xrange, (-Inf, Inf));
    Dxf = dws("x", "f", (N[1],1), xrange, (-Inf, Inf));

    if pol == "TM"
        A = omega^2*mu0*T_eps + Dxf*Dxb;
    elseif pol == "TE"
        A = omega^2*mu0*T_eps + T_eps*Dxf*T_eps_x.^-1*Dxb;
    else
        error("Invalid polarization specified!");
    end

    (betas2, vectors) = eigs(A, nev=Neigs, sigma=beta_est^2);
    betas = sqrt.(betas2+0im);
    return (betas, vectors)
end

function solve_eigen_Bloch_TM(omega_est, Kx, Neigs, xrange, yrange, eps_r, Npml)
    N = size(eps_r);
    M = prod(N);

    T_eps_z = spdiagm(epsilon0*eps_r[:]);

    Hx = zeros(Complex128, Neigs, N[1], N[2]);
    Hy = zeros(Complex128, Neigs, N[1], N[2]);
    Ez = zeros(Complex128, Neigs, N[1], N[2]);

    (Sxf, Sxb, Syf, Syb) = S_create(omega_est, N, Npml, xrange, yrange);

    # Construct derivates
    Dxb = Sxb*dws("x", "b", N, xrange, yrange);
    Dxf = Sxf*dws("x", "f", N, xrange, yrange);
    Dyb = Syb*dws("y", "b", N, xrange, yrange);
    Dyf = Syf*dws("y", "f", N, xrange, yrange);

    A = -mu0^-1*T_eps_z.^-1*(Dxb*Dxf - 1im*Kx*Dxb - 1im*Kx*Dxf - Kx^2*speye(M) + Dyb*Dyf); 

    (omegas_sqr, uz) = eigs(A, nev=Neigs, sigma=omega_est^2);
    omegas = sqrt.(omegas_sqr);

    x_vec = transpose(linspace(xrange(1), xrange(2), N(1))); 
    x_space = repmat(x_vec, [N(2), 1]); 

    for i = 1:Neigs
        ez_temp = uz[:, i].*exp(-1im*Kx*x_space); 
        hx_temp = -1/1im/omega[i]*mu0^-1*Dyf*ez_temp; 
        hy_temp = 1/1im/omega[i]*mu0^-1*Dxf*ez_temp; 
        
        Ez[i, :, :] = reshape(ez_temp, N); 
        Hx[i, :, :] = reshape(hx_temp, N); 
        Hy[i, :, :] = reshape(hy_temp, N); 
    end

    return (omega, Ez, Hx, Hy)
end