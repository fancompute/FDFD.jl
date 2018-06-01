export ϵ₀, μ₀, c₀, η₀
export Polarization, TM, TE
export Direction, x̂, ŷ, ẑ
export Point

const ϵ₀ = 8.85418782e-12;
const μ₀ = 1.25663706e-6;
const c₀ = sqrt(1/ϵ₀/μ₀);
const η₀ = sqrt(μ₀/ϵ₀);

const Float  =  Float64;
const Tuple2 =  NTuple{2};

const XX = 1
const YY = 2
const ZZ = 3

@enum Direction x̂=1 ŷ=2 ẑ=3
@enum DerivativeDirection Forward=1 Backward=2
@enum Polarization TM=1 TE=2
@enum FDFDMatSymmetry CSym=1 CNSym=2

struct Point <: FieldVector{2, Float64}
    x::Float64
    y::Float64
end
