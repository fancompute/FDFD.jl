export solve

function solve(d::Device)
    ω = d.ω[1];

    Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]);

    Hx = zeros(Complex128, size(d.grid));
    Hy = zeros(Complex128, size(d.grid));
    Ez = zeros(Complex128, size(d.grid));

    (Sxf, Sxb, Syf, Syb) = S_create(d.grid, ω);

    # Construct derivates
    δxb = Sxb*δ(DirectionX, Backward, d.grid);
    δxf = Sxf*δ(DirectionX, Forward,  d.grid);
    δyb = Syb*δ(DirectionY, Backward, d.grid);
    δyf = Syf*δ(DirectionY, Forward,  d.grid);

    # Construct system matrix
    A = δxf*μ₀^-1*δxb + δyf*μ₀^-1*δyb + ω^2*Tϵ;
    b = 1im*ω*d.src[:];

    ez = dolinearsolve(A, b, matrixtype=Pardiso.COMPLEX_SYM)

    hx = -1/1im/ω/μ₀*δyb*ez;
    hy = 1/1im/ω/μ₀*δxb*ez;

    Hx = reshape(hx, size(d.grid));
    Hy = reshape(hy, size(d.grid));
    Ez = reshape(ez, size(d.grid));

    return (Ez, Hx, Hy)
end