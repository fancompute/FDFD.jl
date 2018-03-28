export solve_eigen_1D

function solve_eigen_1D(geom::Geometry1D, polarization, ω, estimatedβ, neigenvalues)
    eps_x = grid_average(geom.epsr, "x");
    T_eps = spdiagm(ϵ₀*geom.epsr[:]);
    T_eps_x = spdiagm(ϵ₀*eps_x[:]);

    δxb = δ("x", "b", geom);
    δxf = δ("x", "f", geom);

    if polarization == "TM"
        A = ω^2*μ₀*T_eps + δxf*δxb;
    elseif polarization == "TE"
        A = ω^2*μ₀*T_eps + T_eps*δxf*T_eps_x.^-1*δxb;
    else
        error("Invalid polarization specified!");
    end

    (β², vectors) = eigs(A, nev=neigenvalues, sigma=estimatedβ^2);
    β = sqrt.(β²);
    return (β, vectors)
end
