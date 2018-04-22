export solve_eigen_1D

function solve_eigen_1D(geom::Geometry1D, polarization, ω, estimatedβ, neigenvalues)
    Tϵ = spdiagm(ϵ₀*geom.ϵᵣ[:]);
    Tϵxinv = spdiagm(ϵ₀*grid_average(geom.ϵᵣ, "x")[:].^-1);

    δxb = δ("x", "b", geom);
    δxf = δ("x", "f", geom);

    if polarization == "TM"
        A = ω^2*μ₀*Tϵ + δxf*δxb;
    elseif polarization == "TE"
        A = ω^2*μ₀*Tϵ + Tϵ*δxf*Tϵxinv*δxb;
    else
        error("Invalid polarization specified!");
    end

    (β², vectors) = eigs(A, nev=neigenvalues, sigma=estimatedβ^2);
    β = sqrt.(β²);
    return (β, vectors)
end
