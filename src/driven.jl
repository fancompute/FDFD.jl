export solve_driven_TM

function solve_driven_TM(geom::Geometry2D, ω)
    Tϵ = spdiagm(ϵ₀*geom.ϵᵣ[:]);

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
    A = δxf*μ₀^-1*δxb + δyf*μ₀^-1*δyb + ω^2*Tϵ;
    b = 1im*ω*geom.src[:];

    ez = dolinearsolve(A, b, matrixtype=Pardiso.COMPLEX_SYM)

    hx = -1/1im/ω/μ₀*δyb*ez;
    hy = 1/1im/ω/μ₀*δxb*ez;

    Hx = reshape(hx, geom.N);
    Hy = reshape(hy, geom.N);
    Ez = reshape(ez, geom.N);

    return (Ez, Hx, Hy)
end