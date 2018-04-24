export solve

function solve(d::AbstractDevice, pol::Polarization, estimatedβ::Complex, neigenvalues::Int)
    ω = d.ω[1];
    if ndims(d.grid) == 1
        Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]);
        Tϵxinv = spdiagm(ϵ₀*grid_average(d.ϵᵣ, DirectionX)[:].^-1);

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

    (β², vectors) = eigs(A, nev=neigenvalues, sigma=estimatedβ^2);
    β = sqrt.(β²);
    return (β, vectors)
end