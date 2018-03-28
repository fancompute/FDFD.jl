export solve_driven_TM

function solve_driven_TM(geom::Geometry2D, ω)
    T_eps_z = spdiagm(epsilon0*geom.epsr[:]);

    Hx = zeros(Complex128, geom.N);
    Hy = zeros(Complex128, geom.N);
    Ez = zeros(Complex128, geom.N);

    (Sxf, Sxb, Syf, Syb) = S_create(ω, geom.N, geom.Npml, geom.xrange, geom.yrange);

    # Construct derivates
    δxb = Sxb*δ("x", "b", geom);
    δxf = Sxf*δ("x", "f", geom);
    δyb = Syb*δ("y", "b", geom);
    δyf = Syf*δ("y", "f", geom);

    # Construct system matrix
    A = δxf*μ₀^-1*δxb + δyf*μ₀^-1*δyb + ω^2*T_eps_z;
    b = 1im*ω*geom.src[:];

    if solver_pardiso
        ez = solve(handle_ps, A, b);
    else
        ez = lufact(A)\b;
    end

    hx = -1/1im/ω/μ₀*δyb*ez;
    hy = 1/1im/ω/μ₀*δxb*ez;

    Hx = reshape(hx, geom.N);
    Hy = reshape(hy, geom.N);
    Ez = reshape(ez, geom.N);

    return (Ez, Hx, Hy)
end