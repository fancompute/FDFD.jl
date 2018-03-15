export solve_driven_TM

function solve_driven_TM(geom, omega)
    T_eps_z = spdiagm(epsilon0*geom.epsr[:]);

    Hx = zeros(Complex128, geom.N);
    Hy = zeros(Complex128, geom.N);
    Ez = zeros(Complex128, geom.N);

    (Sxf, Sxb, Syf, Syb) = S_create(omega, geom.N, geom.Npml, geom.xrange, geom.yrange);

    # Construct derivates
    Dxb = Sxb*dws("x", "b", geom);
    Dxf = Sxf*dws("x", "f", geom);
    Dyb = Syb*dws("y", "b", geom);
    Dyf = Syf*dws("y", "f", geom);

    # Construct system matrix
    A = Dxf*mu0^-1*Dxb + Dyf*mu0^-1*Dyb + omega^2*T_eps_z;
    b = 1im*omega*geom.src[:];

    if solver_pardiso
        ez = solve(handle_ps, A, b);
    else
        ez = lufact(A)\b;
    end

    hx = -1/1im/omega/mu0*Dyb*ez;
    hy = 1/1im/omega/mu0*Dxb*ez;

    Hx = reshape(hx, geom.N);
    Hy = reshape(hy, geom.N);
    Ez = reshape(ez, geom.N);

    return (Ez, Hx, Hy)
end