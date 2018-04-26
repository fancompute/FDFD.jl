export ϵ₀, μ₀, c₀, η₀
export Polarization, TM, TE
export Direction, DirectionX, DirectionY

const ϵ₀ = 8.85418782e-12;
const μ₀ = 1.25663706e-6;
const c₀ = sqrt(1/ϵ₀/μ₀);
const η₀ = sqrt(μ₀/ϵ₀);

const Float    =  Float64;
const Tuple2   =  NTuple{2};

@enum Direction 			DirectionX=1 	DirectionY=2
@enum DerivativeDirection 	Forward=1 		Backward=2
@enum Polarization 			TM=1 			TE=2