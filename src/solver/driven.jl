export solve

"solve(d::Device, pol::Polarization)"
function solve(d::Device, pol::Polarization=TM)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);

    Nω = length(d.ω);
    pol == TM && ( fields = Array{FieldTM}(Nω) );
    pol == TE && ( fields = Array{FieldTE}(Nω) );

    for i in eachindex(d.ω)
        print_info("======= Frequency: $i/$Nω =======");
        ω = d.ω[i];

        length(d.modes) > 0 && ( d.src=zeros(Complex, size(d.grid)) ) # Only reset if we are using modes
        for mode in d.modes
            # TODO: handle the polarization here
            setup_mode!(d, TM, ω, mode.neff, mode.coor, mode.dir, mode.width);
        end

        Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]);
        Tϵxi = spdiagm(1./grid_average(ϵ₀*d.ϵᵣ, x̂)[:]);
        Tϵyi = spdiagm(1./grid_average(ϵ₀*d.ϵᵣ, ŷ)[:]);

        (Sxf, Sxb, Syf, Syb) = S_create(d.grid, ω);

        # Construct derivates
        δxb = Sxb*δ(x̂, Backward, d.grid);
        δxf = Sxf*δ(x̂, Forward,  d.grid);
        δyb = Syb*δ(ŷ, Backward, d.grid);
        δyf = Syf*δ(ŷ, Forward,  d.grid);

        if pol == TM
            # Construct system matrix
            A = δxf*μ₀^-1*δxb + δyf*μ₀^-1*δyb + ω^2*Tϵ;
            b = 1im*ω*d.src[:];

            ez = dolinearsolve(A, b, CSym)

            hx = -1/1im/ω/μ₀*δyb*ez;
            hy = 1/1im/ω/μ₀*δxb*ez;
            fields[i] = FieldTM(d.grid, ω, ez, hx, hy);
        elseif pol == TE
            # Construct system matrix
            A = δxf*Tϵxi*δxb + δyf*Tϵyi*δyb + ω^2*μ₀*speye(length(d.grid));
            b = 1im*ω*d.src[:];

            hz = dolinearsolve(A, b, CSym)

            ex = 1/1im/ω*Tϵyi*δyb*hz;
            ey = 1/1im/ω*Tϵxi*(-δxb*hz);

            fields[i] = FieldTE(d.grid, ω, hz, ex, ey);
        end
    end

    Nω == 1 && return fields[1]
    return fields
end
