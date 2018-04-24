export Grid
export dx, dy, dh, xc, yc, xe, ye, coord2ind

# Base Grid struc
struct Grid{K}
    L::SVector{K,Float} # Length of dimensions
    N::SVector{K,Int} # Number of cells
    Npml::SVector{K,Int} # Number of PML cells on each end for each dimension
    bounds::Tuple2{SVector{K,Float}} # Boundary coordinates, in bounds[edge][dim]
end


Base.ndims(::Grid{K}) where {K} = K
Base.length(g::Grid) = prod(g.N)
Base.size(g::Grid{2}) = (g.N[1], g.N[2])
Base.size(g::Grid, i::Int) = g.N[i]
Base.size(g::Grid{1}) = (g.N[1])


# 2D Grid 
function Grid(dh::Float, Npml::AbstractArray{Int,1}, xrange::AbstractArray{Float,1}, yrange::AbstractArray{Float,1})
    L = SVector{2}([xrange[2]-xrange[1], yrange[2]-yrange[1]]);
    N = SVector{2}(Int.(round.(L/dh)));
    bounds1 = SVector{2}([xrange[1], yrange[1]]);
    bounds2 = SVector{2}([xrange[2], yrange[2]]);
    return Grid(L, N, SVector{2}(Npml), (bounds1, bounds2));
end


# 1D Grid 
function Grid(dh::Float, Npml::AbstractArray{Int,1}, xrange::AbstractArray{Float,1})
    L = SVector{1}(xrange[2]-xrange[1]);
    N = SVector{1}(Int.(round.(L/dh)));
    bounds1 = SVector{1}(xrange[1]);
    bounds2 = SVector{1}(xrange[2]);
    return Grid(L, N, SVector{1}(Npml), (bounds1, bounds2));
end

function Grid(N::Int, xrange::AbstractArray{Float,1})
    L = SVector{1}(xrange[2]-xrange[1]);
    N = SVector{1}(N);
    bounds1 = SVector{1}(xrange[1]);
    bounds2 = SVector{1}(xrange[2]);
    return Grid(L, SVector{1}(N), SVector{1}(0), (bounds1, bounds2));
end


# Coordinate functions
function dh(g::Grid, w::Direction)
	w == DirectionX && return dx(g)
    w == DirectionY && return dy(g)
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

function ye(g::Grid)
	return g.bounds[1][2]+dy(g)*(0:1:g.N[2])
end

xe(g::Grid, i::CartesianIndex{2}) = xe(g)[i.I[1]];
ye(g::Grid, i::CartesianIndex{2}) = ye(g)[i.I[2]];
xe(g::Grid, i::Int64) = xe(g)[ ind2sub(size(g),i)[1] ];
ye(g::Grid, i::Int64) = ye(g)[ ind2sub(size(g),i)[2] ];

function coord2ind(g::Grid{K}, xy::AbstractArray) where {K}
    indx = Int(round((xy[1]-g.bounds[1][1])/g.L[1]*size(g,1))+1);
    if K == 1
    	return indx
    else
    	indy = Int(round((xy[2]-g.bounds[1][2])/g.L[2]*size(g,2))+1);
	end
    return (indx, indy)
end