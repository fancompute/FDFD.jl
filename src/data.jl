export Result, Field, FieldTM, FieldTE, FieldAll
export AxisX, AxisY, AxisComponent

const AxisX = Axis{:x}
const AxisY = Axis{:y}
const AxisComponent = Axis{:component, Array{Symbol,1}}

abstract type Result{T,N} <: AbstractArray{T,N} end
abstract type Field{T,N} <: Result{T,N} end

Base.size(A::Result) = size(A.data)
Base.IndexStyle(A::Result) = IndexStyle(A.data)

AxisArrays.axisdim(A::Result, ax) = axisdim(A.data, ax)
AxisArrays.axes(A::Result, i...) = AxisArrays.axes(A.data, i...)
AxisArrays.axisnames(A::Result) = axisnames(A.data)
AxisArrays.axisvalues(A::Result) = axisvalues(A.data)

# function Base.similar(x::Result, ::Type{S}, dims::NTuple{N,Int}) where {S,N}
#     if size(x) == dims
#         data = similar(AxisArray(x),S)
#         similar_helper(x,data,x.params)
#     else
#         similar(AxisArray(x),S,dims)
#     end
# end
#
# Base.similar(x::Result, ax1::Axis, axs::Axis...) = similar(x, eltype(x), ax1, axs...)
#
# function Base.similar(x::Result, ::Type{S}, ax1::Axis, axs::Axis...) where S
#     if ndims(x) == 1+length(axs)
#         data = similar(AxisArray(x), ax1, axs...)
#         similar_helper(x, data, x.params)
#     else
#         similar(AxisArray(x), S, ax1, axs...)
#     end
# end

@inline @Base.propagate_inbounds Base.getindex(A::Result, i...) =
  getindex(A.data, i...)
@inline @Base.propagate_inbounds Base.getindex(A::Result, i::Int...) =
  getindex(A.data, i...)
@inline @Base.propagate_inbounds Base.setindex!(A::Result, val, i...) =
  setindex!(A.data, val, i...)
@inline @Base.propagate_inbounds Base.setindex!(A::Result{T}, val::T, i::Int...) where T =
  setindex!(A.data, val, i...)

# ============================================================================ #

struct FieldTM <: Field{Complex,3}
    grid::Grid{2}
    ω::Complex
    data::AxisArray
end

function FieldTM(grid::Grid{2}, ω::Number, data::Array{<:Complex,3})
    return FieldTM(grid, Complex(ω), AxisArray(data, AxisX(xc(grid)), AxisY(yc(grid)), AxisComponent([:Ez, :Hx, :Hy])))
end

function FieldTM(grid::Grid{2}, ω::Number, Ez::Array{<:Complex,1}, Hx::Array{<:Complex,1}, Hy::Array{<:Complex,1})
  	sz = size(grid);
    return FieldTM(grid, ω, cat(reshape(Ez, sz), reshape(Hx, sz), reshape(Hy, sz), dims=3))
end

struct FieldTE <: Field{Complex,3}
    grid::Grid{2}
    ω::Complex
    data::AxisArray
end

function FieldTE(grid::Grid{2}, ω::Number, data::Array{<:Complex,3})
    return FieldTE(grid, Complex(ω), AxisArray(data, AxisX(xc(grid)), AxisY(yc(grid)), AxisComponent([:Hz, :Ex, :Ey])))
end

function FieldTE(grid::Grid{2}, ω::Number, Hz::Array{<:Complex,1}, Ex::Array{<:Complex,1}, Ey::Array{<:Complex,1})
  	sz = size(grid);
    return FieldTE(grid, ω, cat(reshape(Hz, sz), reshape(Ex, sz), reshape(Ey, sz), dims=3))
end

struct FieldAll <: Field{Complex,3}
    grid::Grid{2}
    ω::Complex
    data::AxisArray
end

function FieldAll(grid::Grid{2}, ω::Number, data::Array{<:Complex,3})
    return FieldAll(grid, Complex(ω), AxisArray(data, AxisX(xc(grid)), AxisY(yc(grid)), AxisComponent([:Ex, :Ey, :Ez, :Hx, :Hy, :Hz])))
end

function FieldAll(grid::Grid{2}, ω::Number, Ex::Array{<:Complex,1}, Ey::Array{<:Complex,1}, Ez::Array{<:Complex,1}, Hx::Array{<:Complex,1}, Hy::Array{<:Complex,1}, Hz::Array{<:Complex,1})
  	sz = size(grid);
    return FieldAll(grid, ω, cat(3, reshape(Ex, sz), reshape(Ey, sz), reshape(Ez, sz), reshape(Hx, sz), reshape(Hy, sz), reshape(Hz, sz)))
end

# ============================================================================ #

struct Flux <: Field{Float,3}
    grid::Grid{2}
    ω::Complex
    data::AxisArray
end

function Flux(grid::Grid{2}, ω::Number, data::Array{Float,3})
    return Flux(grid, Complex(ω), AxisArray(data, AxisX(xc(grid)), AxisY(yc(grid)), AxisComponent([:Sx, :Sy])))
end

function Flux(grid::Grid{2}, ω::Number, Sx::Array{Float,2}, Sy::Array{Float,2})
    return Flux(grid, ω, cat(3, Sx, Sy))
end

# poynt = poynting(field_dipole);
# ax = plot_field(field_dipole, funcz=real)
# xinterval = -2..2
# yinterval = 1..3
# sub_poynt = poynt[xinterval,yinterval,:]
# x_vals = AxisArrays.axisvalues(sub_poynt)[1]
# y_vals = AxisArrays.axisvalues(sub_poynt)[2]
# skip = 9
# quiver(x_vals[1:skip:end], y_vals[1:skip:end], sub_poynt[1:skip:end,1:skip:end,:Sx]', sub_poynt[1:skip:end,1:skip:end,:Sy]')
