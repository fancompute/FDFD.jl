export AbstractDevice, Device
export setup_ϵᵣ!, setup_src!

abstract type AbstractDevice end

Base.size(d::AbstractDevice) = size(d.grid);
Base.size(d::AbstractDevice, i::Int) = size(d.grid, i);
Base.length(d::AbstractDevice) = length(d.grid);

mutable struct Device{D} <: AbstractDevice
	grid::Grid{D}
	ϵᵣ::Array{Complex}
    src::Array{Complex}
    ω::AbstractVector{Float}
end

function Device(grid::Grid, ω::AbstractVector{Float})
	ϵᵣ = ones(Complex, size(grid));
    src = zeros(Complex, size(grid));
	return Device(grid, ϵᵣ, src, ω)
end

Device(grid::Grid, ω::Float) = Device(grid, [ω]);

normalize_parameters(g::Grid) = (ϵ₀*g.L₀, μ₀*g.L₀, c₀/g.L₀);
normalize_parameters(d::AbstractDevice) = normalize_parameters(d.grid);

function _compose_shapes!(pixels::AbstractArray, grid::Grid, shapes::AbstractVector{<:Shape})
    kd = KDTree(shapes);
    for i in eachindex(pixels)
        x = xc(grid, i);
        y = yc(grid, i);
        x1 = xe(grid, i);
        y1 = ye(grid, i);
        x2 = xe(grid, i+1);
        y2 = ye(grid, i+1);
        shape = findin([x,y], kd);
        if ~isnull(shape)
            r₀, nout = surfpt_nearby( [x,y], get(shape));
            frac = volfrac((SVector(x1, y1), SVector(x2, y2)), nout, r₀);
            data = get(shape).data;
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

setup_ϵᵣ!(d::AbstractDevice, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.ϵᵣ, d.grid, shapes)
setup_ϵᵣ!(d::AbstractDevice, region, value) = _mask_values!(d.ϵᵣ, d.grid, region, value)

setup_src!(d::AbstractDevice, region, value) = _mask_values!(d.src, d.grid, region, value)

function setup_src!(d::AbstractDevice, xy::AbstractArray) 
    (indx, indy) = coord2ind(d.grid, xy);
    d.src[indx, indy] = 1im;
end

function setup_src!(d::AbstractDevice, pol::Polarization, neff::Number, srcxy::AbstractArray, srcnormal::Direction, srcpoints::Int)
    (indx, indy) = coord2ind(d.grid, srcxy);
    M = Int64(round((srcpoints-1)/2));
    srcpoints = 2*M+1;

    if srcnormal == DirectionX
        indx = indx;
        indy = indy+(-M:M);
        dh = dy(d.grid);
    elseif srcnormal == DirectionY
        indx = indx+(-M:M);
        indy = indy;
        dh = dx(d.grid);
    end

    g1D   = Grid(srcpoints, [0 srcpoints*dh], L₀=d.grid.L₀);
    dev1D = Device(g1D, d.ω);
    dev1D.ϵᵣ = d.ϵᵣ[indx, indy];
    (β, vector) = solve(dev1D, pol, neff, 1);
    d.src[indx, indy] = abs.(vector);
end
