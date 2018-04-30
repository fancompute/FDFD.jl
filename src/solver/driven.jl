export solve

"    solve(d::Device)"
function solve(d::Device)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);

    Nω = length(d.ω);
    fields = Array{FieldTM}(Nω);

    for i in eachindex(d.ω)
        ω = d.ω[i];

        length(d.modes) > 0 && ( d.src=zeros(Complex, size(d.grid)) ) # Only reset if we are using modes
        for mode in d.modes
            # TODO: handle the polarization here
            setup_mode!(d, TM, ω, mode.neff, mode.coor, mode.dir, mode.width);
        end

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
        fields[i] = FieldTM(d.grid, ez, hx, hy);
    end
    
    Nω == 1 && return fields[1]
    return fields
end