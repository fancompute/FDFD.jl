export SolvedResult, Field, FieldTM, FieldTE, Flux2D, FieldSlice

abstract type SolvedResult end
abstract type Field <: SolvedResult end

mutable struct Flux2D <: SolvedResult
    grid::Grid{2}
    Sx::Array{Float,2}
    Sy::Array{Float,2}
end

mutable struct FieldTM <: Field
    grid::Grid{2}
    Ez::Array{Complex,2}
    Hx::Array{Complex,2}
    Hy::Array{Complex,2}
end

function FieldTM(grid::Grid{2}, Ez::Array{<:Complex,1}, Hx::Array{<:Complex,1}, Hy::Array{<:Complex,1})
	sz = size(grid);
	return FieldTM(grid, reshape(Ez, sz), reshape(Hx, sz), reshape(Hy, sz))
end

mutable struct FieldTE <: Field
    grid::Grid{2}
    Hz::Array{Complex,2}
    Ex::Array{Complex,2}
    Ey::Array{Complex,2}
end

mutable struct FieldSlice <: Field
    grid::Grid{1}
    val::Array{Complex,1}
end