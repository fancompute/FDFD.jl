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

    (beta_sqr, output_vectors) = eigs(A, nev=Neigs, sigma=beta_est^2);
    beta = sqrt.(beta_sqr);
    return (beta, output_vectors)
end