export eigenmode, eigenfrequency

"    eigenmode(d::AbstractDevice, pol::Polarization, neff::Number, neigenvalues::Int)"
function eigenmode(d::AbstractDevice, pol::Polarization, neff::Number, neigenvalues::Int)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);
    ω = d.ω[1];
    if ndims(d.grid) == 1
        Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]);
        Tϵxinv = spdiagm((ϵ₀*grid_average(d.ϵᵣ, x̂)[:]).^-1);

        δxb = δ(x̂, Backward, d.grid);
        δxf = δ(x̂, Forward,  d.grid);

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

"    eigenfrequency(d::AbstractDevice, pol::Polarization, neigenvalues::Int; which::Symbol=:LM)"
function eigenfrequency(d::AbstractDevice, pol::Polarization, neigenvalues::Int; which::Symbol=:LM)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);
    ω₀ = d.ω[1];

    (Sxf, Sxb, Syf, Syb) = S_create(d.grid, ω₀);

    δxb = Sxb*δ(x̂, Backward, d.grid);
    δxf = Sxf*δ(x̂, Forward,  d.grid);
    δyb = Syb*δ(ŷ, Backward, d.grid);
    δyf = Syf*δ(ŷ, Forward,  d.grid);

    if pol == TM
        fields = Array{FieldTM}(neigenvalues);

        Tϵᵣ⁻¹ = spdiagm(1./d.ϵᵣ[:]);
        A = Tϵᵣ⁻¹*δxf*δxb + Tϵᵣ⁻¹*δyf*δyb;

        (ω²μϵ, ez) = eigs(A, nev=neigenvalues, sigma=-ω₀^2*μ₀*ϵ₀, which=which);
        ω = sqrt.(-ω²μϵ/μ₀/ϵ₀);

        for i = 1:neigenvalues
            hxi = -1/1im/ω[i]/μ₀*δyf*ez[:,i];
            hyi = 1/1im/ω[i]/μ₀*δxf*ez[:,i];
            fields[i] = FieldTM(d.grid, ω[i], ez[:,i], hxi, hyi)
        end

        return (ω, fields)
    end
    if pol == TE
        fields = Array{FieldTE}(neigenvalues);

        Tϵx⁻¹ = spdiagm(1./grid_average(ϵ₀*d.ϵᵣ, x̂)[:]);
        Tϵy⁻¹ = spdiagm(1./grid_average(ϵ₀*d.ϵᵣ, ŷ)[:]);
        A = δxf*Tϵx⁻¹*δxb + δyf*Tϵy⁻¹*δyb;

        (ω²μ, hz) = eigs(A, nev=neigenvalues, sigma=-ω₀^2*μ₀, which=which);
        ω = sqrt.(-ω²μ/μ₀);

        for i = 1:neigenvalues
            exi = 1/1im/ω[i]*Tϵx⁻¹*δyb*hz[:,i];
            eyi = 1/1im/ω[i]*Tϵy⁻¹*(-δxb*hz[:,i]);
            fields[i] = FieldTE(d.grid, ω[i], hz[:,i], exi, eyi)
        end

        return (ω, fields)
    end
end
