export solve

"    solve(d::AbstractDevice, pol::Polarization, neff::Number, neigenvalues::Int)"
function solve(d::AbstractDevice, pol::Polarization, neff::Number, neigenvalues::Int)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);
    ω = d.ω[1];
    if ndims(d.grid) == 1
        Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]);
        Tϵxinv = spdiagm((ϵ₀*grid_average(d.ϵᵣ, DirectionX)[:]).^-1);

        δxb = δ(DirectionX, Backward, d.grid);
        δxf = δ(DirectionX, Forward,  d.grid);

        if pol == TM
            A = ω^2*μ₀*Tϵ + δxf*δxb;
        end
        if pol == TE
            A = ω^2*μ₀*Tϵ + Tϵ*δxf*Tϵxinv*δxb;
        end
    else
        error("solve() currently only solves 1D systems");
    end

    estimatedβ = ω/c₀*neff;

    (β², vectors) = eigs(A, nev=neigenvalues, sigma=estimatedβ^2);

    β = sqrt.(β²);
    return (β, vectors)
end