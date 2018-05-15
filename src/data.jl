export SolvedResult, Field, FieldTM, FieldTE, Flux2D, FieldSlice

abstract type SolvedResult end
abstract type Field <: SolvedResult end

struct Flux2D <: SolvedResult
    grid::Grid{2}
    Sx::Array{Float,2}
    Sy::Array{Float,2}
end

struct FieldTM <: Field
    grid::Grid{2}
    Ez::Array{Complex,2}
    Hx::Array{Complex,2}
    Hy::Array{Complex,2}
end

function FieldTM(grid::Grid{2}, Ez::Array{<:Complex,1}, Hx::Array{<:Complex,1}, Hy::Array{<:Complex,1})
	sz = size(grid);
	return FieldTM(grid, reshape(Ez, sz), reshape(Hx, sz), reshape(Hy, sz))
end

struct FieldTE <: Field
    grid::Grid{2}
    Hz::Array{Complex,2}
    Ex::Array{Complex,2}
    Ey::Array{Complex,2}
end

function FieldTE(grid::Grid{2}, Hz::Array{<:Complex,1}, Ex::Array{<:Complex,1}, Ey::Array{<:Complex,1})
    sz = size(grid);
    return FieldTE(grid, reshape(Hz, sz), reshape(Ex, sz), reshape(Ey, sz))
end

struct FieldSlice <: Field
    grid::Grid{1}
    val::Array{Complex,1}
end