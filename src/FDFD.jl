module FDFD

using GeometryPrimitives, StaticArrays, AxisArrays, IntervalSets
using Logging, SparseArrays, LinearAlgebra

include("./types.jl");
include("./grid.jl");
include("./data.jl");
include("./device.jl");
include("./pml.jl");
include("./flux.jl");
include("./solver/solver.jl");
include("./solver/driven.jl");
include("./solver/eigen.jl");
include("./solver/modulation.jl");
include("./solver/nonlinear.jl");

export unwrap, set_log_level!

global log_level = 1

set_log_level!(level::Int) = global log_level = level

iscallable(f) = !isempty(methods(f))

function unwrap(v, inplace=false)
  # currently assuming an array
  unwrapped = inplace ? v : copy(v)
  for i in 2:length(v)
    while unwrapped[i] - unwrapped[i-1] >= pi
      unwrapped[i] -= 2pi
    end
    while unwrapped[i] - unwrapped[i-1] <= -pi
      unwrapped[i] += 2pi
    end
  end
  return unwrapped
end

end
