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
end

"    χ3Device(grid::Grid, ω::AbstractVector{Float})"
function χ3Device(grid::Grid, ω::AbstractVector{Float})
    ϵᵣ = ones(Complex, size(grid));
    χ = zeros(Float, size(grid));
    src = zeros(Complex, size(grid));
    return χ3Device(grid, ϵᵣ, χ, src, ω)
end

"    setup_χ!(d::χ3Device, shapes::AbstractVector{<:Shape})"
setup_χ!(d::χ3Device, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.χ, d.grid, shapes)

"    setup_χ!(d::χ3Device, region, value)"
setup_χ!(d::χ3Device, region, value) = _mask_values!(d.χ, d.grid, region, value)

"    solve(d::χ3Device, which_method::IterativeMethod)"
function solve(d::χ3Device, which_method::IterativeMethod)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);
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

    A = δxf*μ₀^-1*δxb + δyf*μ₀^-1*δyb + ω^2*Tϵ;
    b = 1im*ω*d.src[:];
    print_info("Solving linear system");
    ez = dolinearsolve(A, b, CNSym);

    coeff = ω^2*ϵ₀*3*d.χ[:]/d.grid.L₀;

    if which_method == IterativeMethodB
        print_info("Starting nonlinear iteration using Born");
        (ez, err) = _doborn(ez, A, b, coeff);
    end
    if which_method == IterativeMethodGN
        print_info("Starting nonlinear iteration using Gauss-Newton");
        (ez, err) = _donewton(ez, A, b, coeff);
    end

    hx = -1/1im/ω/μ₀*δyb*ez;
    hy = 1/1im/ω/μ₀*δxb*ez;
    
    return (FieldTM(d.grid, ez, hx, hy), err)
end

function _doborn(ez, A, b, coeff; tol = 1e-12, maxiterations = 50)
    i = 1;
    err = [1.0];
    while err[end] > tol && i <= maxiterations
        print_info(@sprintf("iteration number: %d", i));
        ez_new = dolinearsolve(A + spdiagm(coeff.*ez.*conj.(ez)), b, CNSym);
        append!(err, norm(ez_new - ez)/norm(ez));
        
        ez = ez_new;
        i += 1;

        print_info(@sprintf("step error: %e", err[end]));
    end
    return (ez, err)
end

function _donewton(ez, A, b, coeff; tol = 1e-12, maxiterations = 50)
    i = 1;
    err = [1.0];
    M = length(ez);

    ez = [ez; conj.(ez)];
    while err[end] > tol && i <= maxiterations
        print_info(@sprintf("iteration number: %d", i));
        F = (A + spdiagm(coeff.*ez[1:M].*conj.(ez[1:M])))*ez[1:M] - b;
        J1  = A + spdiagm(2*coeff.*conj.(ez[1:M]).*ez[1:M]);
        J2 = spdiagm(coeff.*ez[1:M].*ez[1:M]);
        J = [J1 J2; conj.(J2) conj.(J1)];
        δez = dolinearsolve(J, [F; conj.(F)], CNSym);
        normez = norm(ez);
        ez = ez - δez;
        append!(err, norm(δez)/normez);
        i += 1;
        print_info(@sprintf("step error: %e", err[end]));
    end
    ez = ez[1:M];
    return (ez, err)
end
