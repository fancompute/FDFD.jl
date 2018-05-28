export Grid
export dx, dy, dh, xc, yc, xe, ye
export coord2ind, x2ind, y2ind

const DEFAULT_L₀ = 1e-6;

struct Grid{D}
    L::SVector{D,Float} # Length of dimensions
    L₀::Float # Length unit
    N::SVector{D,Int64} # Number of cells
    Npml::SVector{D,Int64} # Number of PML cells on each end for each dimension
    bounds::Tuple2{SVector{D,Float}} # Boundary coordinates, accessed as bounds[edge][dim]
end

Base.ndims(::Grid{D}) where {D} = D
Base.length(g::Grid) = prod(g.N)
Base.size(g::Grid{2}) = (g.N[1], g.N[2])
Base.size(g::Grid{1}) = (g.N[1])
Base.size(g::Grid{D}, i::Integer) where {D} = i > D ? 1 : g.N[i]

"""
    Grid(dh::Number, Npml::AbstractArray{<:Integer}, xrange::AbstractArray{<:Real}, yrange::AbstractArray{<:Real}; L₀=DEFAULT_L₀)

Constructs a 2D Grid.
"""
function Grid(dh::Number, Npml::AbstractArray{<:Integer}, xrange::AbstractArray{<:Real}, yrange::AbstractArray{<:Real}; L₀=DEFAULT_L₀)
    L = SVector{2}([ Float(xrange[2]-xrange[1]), Float(yrange[2]-yrange[1])]);
    N = SVector{2}(Int.(round.(L/dh)));
    bounds1 = SVector{2}([Float(xrange[1]), Float(yrange[1])]);
    bounds2 = SVector{2}([Float(xrange[2]), Float(yrange[2])]);
    print_info("Grid size: $N");
    return Grid(L, L₀, N, SVector{2}(Npml), (bounds1, bounds2));
end

"""
    Grid(dh::Number, Npml::Integer, xrange::AbstractArray{<:Real}; L₀=DEFAULT_L₀)

Constructs a 1D Grid.
"""
function Grid(dh::Number, Npml::Integer, xrange::AbstractArray{<:Real}; L₀=DEFAULT_L₀)
    L = SVector{1}(Float(xrange[2]-xrange[1]));
    N = SVector{1}(Int.(round.(L/dh)));
    bounds1 = SVector{1}(Float(xrange[1]));
    bounds2 = SVector{1}(Float(xrange[2]));
    return Grid(L, L₀, N, SVector{1}(Npml), (bounds1, bounds2));
end

"""
    Grid(N::Integer, xrange::AbstractArray{<:Real}; L₀=DEFAULT_L₀)

Constructs a 1D Grid.
"""
function Grid(N::Integer, xrange::AbstractArray{<:Real}; L₀=DEFAULT_L₀)
    L = SVector{1}(Float(xrange[2]-xrange[1]));
    N = SVector{1}(N);
    bounds1 = SVector{1}(Float(xrange[1]));
    bounds2 = SVector{1}(Float(xrange[2]));
    return Grid(L, L₀, SVector{1}(N), SVector{1}(0), (bounds1, bounds2));
end


"Returns the cell size for the grid `g` in the direction specified by `w`"
function dh(g::Grid, w::Direction)
	w == x̂ && return dx(g)
    w == ŷ && return dy(g)
end

function dx(g::Grid)
	return (g.bounds[2][1]-g.bounds[1][1])/g.N[1]
end

function dy(g::Grid{2})
	return (g.bounds[2][2]-g.bounds[1][2])/g.N[2]
end

function xc(g::Grid)
	return g.bounds[1][1]+dx(g)*(0.5:1:g.N[1])
end

function yc(g::Grid{2})
	return g.bounds[1][2]+dy(g)*(0.5:1:g.N[2])
end

xc(g::Grid, i::CartesianIndex{2}) = xc(g)[i.I[1]];
yc(g::Grid, i::CartesianIndex{2}) = yc(g)[i.I[2]];
xc(g::Grid, i::Int64) = xc(g)[ ind2sub(size(g),i)[1] ];
yc(g::Grid, i::Int64) = yc(g)[ ind2sub(size(g),i)[2] ];

function xe(g::Grid)
	return g.bounds[1][1]+dx(g)*(0:1:g.N[1])
end

function ye(g::Grid{2})
	return g.bounds[1][2]+dy(g)*(0:1:g.N[2])
end

xe(g::Grid, i::CartesianIndex{2}) = xe(g)[i.I[1]];
ye(g::Grid, i::CartesianIndex{2}) = ye(g)[i.I[2]];
xe(g::Grid, i::Int64) = xe(g)[ ind2sub(size(g),i)[1] ];
ye(g::Grid, i::Int64) = ye(g)[ ind2sub(size(g),i)[2] ];

"    coord2ind(g::Grid{D}, xy::AbstractArray)"
function coord2ind(g::Grid{D}, xy::AbstractArray) where {D}
    indx = x2ind(g, xy[1]);
    if D == 1
    	return indx
    end
    return (indx, y2ind(g, xy[2]))
end

"    x2ind(g::Grid, x::Real)"
function x2ind(g::Grid, x::Real)
    ind = Int(round((x-g.bounds[1][1])/g.L[1]*size(g,1))+1);
    ind < 1 && return 1
    ind > g.N[1] && return g.N[1]
    return ind
end

"    y2ind(g::Grid, y::Real)"
function y2ind(g::Grid, y::Real)
    ind = Int(round((y-g.bounds[1][2])/g.L[2]*size(g,2))+1);
    ind < 1 && return 1
    ind > g.N[2] && return g.N[2]
    return ind
end

"    δ(w::Direction, s::DerivativeDirection, g::Grid)"
function δ(w::Direction, s::DerivativeDirection, g::Grid{K}) where {K}
    if K == 2
        (Nx, Ny) = size(g);
    else
        Nx = size(g, 1);
        Ny = 1;
    end

    if w == x̂
        if s == Forward
            δxf = 1/dx(g)*spdiagm([ones(Nx-1), -ones(Nx), 1], [1, 0, -Nx+1]);
            return kron(speye(Ny), δxf)
        else
            δxb = 1/dx(g)*spdiagm([-ones(Nx-1), ones(Nx), -1], [-1, 0, Nx-1]);
            return kron(speye(Ny), δxb)
        end
    end
    if w == ŷ
        if s == Backward
            δyf = 1/dy(g)*spdiagm([ones(Ny-1), -ones(Ny), 1], [1, 0, -Ny+1]);
            return kron(δyf, speye(Nx))
        else
            δyb = 1/dy(g)*spdiagm([-ones(Ny-1), ones(Ny), -1], [-1, 0, Ny-1]);
            return kron(δyb, speye(Nx))
        end
    end
end


"    grid_average(centerarray::AbstractArray, w::Direction)"
function grid_average(centerarray::AbstractArray, w::Direction)
    ndims(centerarray) == 1 && return (centerarray+circshift(centerarray, (1)))/2
    w == x̂ && return (centerarray+circshift(centerarray, (1, 0)))/2;
    w == ŷ && return (centerarray+circshift(centerarray, (0, 1)))/2;
    return centerarray
end
