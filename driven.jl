export solve_driven_TM

function solve_driven_TM(omega, xrange, yrange, eps_r, Jz, Npml)
    N = size(eps_r);

    T_eps_z = spdiagm(epsilon0*eps_r[:]);

    jz = Jz[:];

    Hx = zeros(Complex128, N);
    Hy = zeros(Complex128, N);
    Ez = zeros(Complex128, N);

    (Sxf, Sxb, Syf, Syb) = S_create(omega, N, Npml, xrange, yrange);

    # Construct derivates
    Dxb = Sxb*dws("x", "b", N, xrange, yrange);
    Dxf = Sxf*dws("x", "f", N, xrange, yrange);
    Dyb = Syb*dws("y", "b", N, xrange, yrange);
    Dyf = Syf*dws("y", "f", N, xrange, yrange);

    # Construct system matrix
    A = Dxf*mu0^-1*Dxb + Dyf*mu0^-1*Dyb + omega^2*T_eps_z;
    b = 1im*omega*jz;

    ez = A\b;

    hx = -1/1im/omega/mu0*Dyb*ez;
    hy = 1/1im/omega/mu0*Dxb*ez;

    Hx = reshape(hx, N);
    Hy = reshape(hy, N);
    Ez = reshape(ez, N);

    return (Ez, Hx, Hy)
end