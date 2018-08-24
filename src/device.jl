export AbstractDevice, Device
export setup_ϵᵣ!, setup_src!
export Mode, add_mode!

struct Mode
    pol::Polarization
    dir::Direction
    neff::Number
    pt::Point
    width::Number
end

abstract type AbstractDevice{D} end

Base.size(d::AbstractDevice) = size(d.grid);
Base.size(d::AbstractDevice, i::Int) = size(d.grid, i);
Base.length(d::AbstractDevice) = length(d.grid);

mutable struct Device{D} <: AbstractDevice{D}
	grid::Grid{D}
	ϵᵣ::Array{Complex}
    src::Array{Complex}
    ω::AbstractVector{Float}
    modes::Array{Mode}
end

"    Device(grid::Grid, ω::AbstractVector{Float})"
function Device(grid::Grid, ω::AbstractVector{Float})
	ϵᵣ = ones(Complex, size(grid));
    src = zeros(Complex, size(grid));
	return Device(grid, ϵᵣ, src, ω, Array{Mode}(undef, 0))
end

"    Device(grid::Grid, ω::Float)"
Device(grid::Grid, ω::Float) = Device(grid, [ω]);

# ============================================================================ #

"    normalize_parameters(g::Grid)"
normalize_parameters(g::Grid) = (ϵ₀*g.L₀, μ₀*g.L₀, c₀/g.L₀);

"    normalize_parameters(d::AbstractDevice)"
normalize_parameters(d::AbstractDevice) = normalize_parameters(d.grid);

# ============================================================================ #

function _compose_shapes!(pixels::AbstractArray, grid::Grid, shapes::AbstractVector{<:Shape})
    kd = KDTree(shapes);
    for i in eachindex(pixels)
        x = xc(grid, i);
        y = yc(grid, i);
        x1 = xe(grid, i);
        y1 = ye(grid, i);
        x2 = xe(grid, i+1);
        y2 = ye(grid, i+1);
        shape = findin([x, y, 0.0], kd);
        if shape.hasvalue
            r₀, nout = surfpt_nearby( [x, y, 0.0], get(shape));
            frac = volfrac((SVector(x1, y1, 0.0), SVector(x2, y2, 0.0)), nout, r₀);
            data = shape.value.data;
            if iscallable(data)
                pixels[i] = data(x, y);
            else
                pixels[i] = data;
            end
        end
    end
end

function _mask_values!(pixels::AbstractArray, grid::Grid, region, value)
    if ndims(grid) == 2
        mask = [region(x, y) for x in xc(grid), y in yc(grid)];
        if iscallable(value)
            value_assigned = [value(x, y) for x in xc(grid), y in yc(grid)];
            pixels[mask]   = value_assigned[mask];
        else
            pixels[mask]   = value;
        end
    elseif ndims(grid) == 1
        mask = [region(x) for x in xc(grid)];
        if iscallable(value)
            value_assigned = [value(x) for x in xc(grid)];
            pixels[mask]   = value_assigned[mask];
        else
            pixels[mask]   = value;
        end
    else
        error("Unkown device dimension!")
    end;
end

"    setup_ϵᵣ!(d::AbstractDevice, shapes::AbstractVector{<:Shape})"
setup_ϵᵣ!(d::AbstractDevice, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.ϵᵣ, d.grid, shapes)

"    setup_ϵᵣ!(d::AbstractDevice, region, value)"
setup_ϵᵣ!(d::AbstractDevice, region, value) = _mask_values!(d.ϵᵣ, d.grid, region, value)

"    setup_src!(d::AbstractDevice, region, value)"
setup_src!(d::AbstractDevice, region, value) = _mask_values!(d.src, d.grid, region, value)

"    setup_src!(d::AbstractDevice, pt::Point)"
function setup_src!(d::AbstractDevice, pt::Point)
    (indx, indy) = coord2ind(d.grid, pt);
    d.src[indx, indy] = 1im;
end

"    setup_src!(d::AbstractDevice, pt::Point, srcnormal::Direction)"
function setup_src!(d::AbstractDevice, pt::Point, srcnormal::Direction)
    (indx, indy) = coord2ind(d.grid, pt);
    if srcnormal == x̂
        d.src[indx, :] = 1im;
    elseif srcnormal == ŷ
        d.src[:, indy] = 1im;
    end
end

# ============================================================================ #

"    add_mode!(d::AbstractDevice, mode::Mode)"
function add_mode!(d::AbstractDevice, mode::Mode)
    append!(d.modes, [mode]);
end

"    setup_mode!(d::AbstractDevice, pol::Polarization, ω::Float, neff::Number, pt::Point, srcnormal::Direction, srcwidth::Number)"
function setup_mode!(d::AbstractDevice, pol::Polarization, ω::Float, neff::Number, pt::Point, srcnormal::Direction, srcwidth::Number)
    (_, vector, indx, indy) = get_modes(d, pol, ω, neff, 1, pt, srcnormal, srcwidth)
    d.src[indx, indy] += normalize(abs.(vector[:]));
end

"    get_modes(d::AbstractDevice, pol::Polarization, ω::Float, neff::Number, nmodes::Int, pt::Point, slicenormal::Direction, slicewidth::Number)"
function get_modes(d::AbstractDevice, pol::Polarization, ω::Float, neff::Number, nmodes::Int, pt::Point, slicenormal::Direction, slicewidth::Number)
    (indx, indy) = coord2ind(d.grid, pt);

    slicenormal == x̂ && (srcpoints = Int(round(slicewidth/dy(d.grid))));
    slicenormal == ŷ && (srcpoints = Int(round(slicewidth/dx(d.grid))));
    srcpoints % 2 == 0 && (srcpoints += 1); # source points is odd

    M = Int((srcpoints-1)/2);
    srcpoints = 2*M+1;

    if slicenormal == x̂
        indx = indx;
        indy = indy+(-M:M);
        dh = dy(d.grid);
    elseif slicenormal == ŷ
        indx = indx+(-M:M);
        indy = indy;
        dh = dx(d.grid);
    end

    g1D   = Grid(srcpoints, [0 srcpoints*dh], L₀=d.grid.L₀);
    dev1D = Device(g1D, ω);
    dev1D.ϵᵣ = d.ϵᵣ[indx, indy];
    (β, vectors) = eigenmode(dev1D, pol, neff, nmodes);
    return (β, vectors, indx, indy)
end
