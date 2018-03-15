
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
