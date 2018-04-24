export Polarization, TM, TE
export Direction, DirectionX, DirectionY

const Float    =  Float64
const Tuple2   =  NTuple{2}

@enum Direction 			DirectionX=1 	DirectionY=2
@enum DerivativeDirection 	Forward=1 		Backward=2
@enum Polarization 			TM=1 			TE=2