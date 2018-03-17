export solve_eigen_1D, solve_eigen_TM_Bloch

function solve_eigen_1D(geom::Geometry1D, pol, omega, beta_est, Neigs)
    eps_x = grid_average(geom.epsr, "x");
    T_eps = spdiagm(epsilon0*geom.epsr[:]);
    T_eps_x = spdiagm(epsilon0*eps_x[:]);

    Dxb = dws("x", "b", geom);
    Dxf = dws("x", "f", geom);

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

function solve_eigen_TM_Bloch(geom::Geometry2D, omega_est, kx, Neigs)
    T_eps_z = spdiagm(epsilon0*geom.epsr[:]);

    Hx = zeros(Complex128, Neigs, geom.N[1], geom.N[2]);
    Hy = zeros(Complex128, Neigs, geom.N[1], geom.N[2]);
    Ez = zeros(Complex128, Neigs, geom.N[1], geom.N[2]);

    (Sxf, Sxb, Syf, Syb) = S_create(omega_est, geom.N, geom.Npml, geom.xrange, geom.yrange);

    Dxb = Sxb*dws("x", "b", geom);
    Dxf = Sxf*dws("x", "f", geom);
    Dyb = Syb*dws("y", "b", geom);
    Dyf = Syf*dws("y", "f", geom);

    A = -mu0^-1*T_eps_z.^-1*(Dxb*Dxf - 1im*kx*Dxb - 1im*kx*Dxf - kx^2*I + Dyb*Dyf); 

    (omegas2, uz) = eigs(A, nev=Neigs, sigma=omega_est^2);
    omegas = sqrt.(omegas2+0im);

    x_space = repmat(xe(geom), [geom.N(2), 1]); 
    for i = 1:Neigs
        ez = uz[:, i].*exp(-1im*kx*x_space); 
        hx = -1/(1im*omegas[i])/mu0*Dyf*ez; 
        hy = 1/(1im*omegas[i])/mu0*Dxf*ez; 
      
        Hx[i, :, :] = reshape(hx, geom.N);
        Hy[i, :, :] = reshape(hy, geom.N);
        Ez[i, :, :] = reshape(ez, geom.N);
    end

    return (omegas, Ez, Hx, Hy)
end