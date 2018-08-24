export AbstractNonlinearDevice, χ3Device
export setup_χ!, solve
export IterativeMethod, IterativeMethodGN, IterativeMethodB

@enum IterativeMethod IterativeMethodGN=1 IterativeMethodB=2

abstract type AbstractNonlinearDevice{D} <: AbstractDevice{D} end

mutable struct χ3Device{D} <: AbstractNonlinearDevice{D}
    grid::Grid{D}
    ϵᵣ::Array{Complex}
    χ::Array{Float}
    src::Array{Complex}
    ω::AbstractVector{Float}
    modes::Array{Mode}
end

"    χ3Device(grid::Grid, ω::AbstractVector{Float})"
function χ3Device(grid::Grid, ω::AbstractVector{Float})
    ϵᵣ = ones(Complex, size(grid));
    χ = zeros(Float, size(grid));
    src = zeros(Complex, size(grid));
    return χ3Device(grid, ϵᵣ, χ, src, ω, Array{Mode}(0))
end

"    χ3Device(grid::Grid, ω::Float)"
χ3Device(grid::Grid, ω::Float) = χ3Device(grid::Grid, [ω])

"    setup_χ!(d::χ3Device, shapes::AbstractVector{<:Shape})"
setup_χ!(d::χ3Device, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.χ, d.grid, shapes)

"    setup_χ!(d::χ3Device, region, value)"
setup_χ!(d::χ3Device, region, value) = _mask_values!(d.χ, d.grid, region, value)

"    solve(d::χ3Device, which_method::IterativeMethod)"
function solve(d::χ3Device, which_method::IterativeMethod)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);

    Nω = length(d.ω);
    fields = Array{FieldTM}(Nω);

    for i in eachindex(d.ω)
        @info "Frequency: $i/$Nω"
        ω = d.ω[i];

        length(d.modes) > 0 && ( d.src=zeros(Complex, size(d.grid)) ) # Only reset if we are using modes
        for mode in d.modes
            # TODO: handle the polarization here
            setup_mode!(d, TM, ω, mode.neff, mode.pt, mode.dir, mode.width);
        end

        Tϵ = Sparse(Diagonal((ϵ₀*d.ϵᵣ[:]));

        Hx = zeros(ComplexF64, size(d.grid));
        Hy = zeros(ComplexF64, size(d.grid));
        Ez = zeros(ComplexF64, size(d.grid));

        (Sxf, Sxb, Syf, Syb) = S_create(d.grid, ω);

        # Construct derivates
        δxb = Sxb*δ(x̂, Backward, d.grid);
        δxf = Sxf*δ(x̂, Forward,  d.grid);
        δyb = Syb*δ(ŷ, Backward, d.grid);
        δyf = Syf*δ(ŷ, Forward,  d.grid);

        A = δxf*μ₀^-1*δxb + δyf*μ₀^-1*δyb + ω^2*Tϵ;
        b = 1im*ω*d.src[:];
        @info "Solving linear system"
        ez = dolinearsolve(A, b, CNSym);

        coeff = ω^2*ϵ₀*3*d.χ[:]/d.grid.L₀;

        if which_method == IterativeMethodB
            @info "Starting nonlinear iteration using Born"
            (ez, err) = _doborn(ez, A, b, coeff);
        end
        if which_method == IterativeMethodGN
            @info "Starting nonlinear iteration using Gauss-Newton"
            (ez, err) = _donewton(ez, A, b, coeff);
        end

        hx = -1/1im/ω/μ₀*δyb*ez;
        hy = 1/1im/ω/μ₀*δxb*ez;

        fields[i] = FieldTM(d.grid, ω, ez, hx, hy)
    end

    Nω == 1 && return fields[1]
    return fields
end

function _doborn(ez, A, b, coeff; tol = 1e-12, maxiterations = 50)
    i = 1;
    err = [1.0];
    while err[end] > tol && i <= maxiterations
        @debug "iteration number: $i"
        ez_new = dolinearsolve(A + Sparse(Diagonal(coeff.*ez.*conj.(ez)), b, CNSym));
        append!(err, norm(ez_new - ez)/norm(ez));

        ez = ez_new;
        i += 1;

        @debug "step error: $(err[end])"
    end
    return (ez, err)
end

function _donewton(ez, A, b, coeff; tol = 1e-12, maxiterations = 50)
    i = 1;
    err = [1.0];
    M = length(ez);

    ez = [ez; conj.(ez)];
    while err[end] > tol && i <= maxiterations
        @debug "iteration number: $i"
        F = (A + Sparse(Diagonal(coeff.*ez[1:M].*conj.(ez[1:M]))))*ez[1:M] - b;
        J1  = A + Sparse(Diagonal(2*coeff.*conj.(ez[1:M]).*ez[1:M]));
        J2 = Sparse(Diagonal(coeff.*ez[1:M].*ez[1:M]));
        J = [J1 J2; conj.(J2) conj.(J1)];
        δez = dolinearsolve(J, [F; conj.(F)], CNSym);
        normez = norm(ez);
        ez = ez - δez;
        append!(err, norm(δez)/normez);
        i += 1;
        @debug "step error: $(err[end])"
    end
    ez = ez[1:M];
    return (ez, err)
end
