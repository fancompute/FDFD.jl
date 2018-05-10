module FDFD

using GeometryPrimitives, StaticArrays

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
include("./plot.jl");

export poynting, flux_surface, unwrap, set_log_level!

global log_level = 1;

function set_log_level!(level::Int)
    global log_level = level
end

iscallable(f) = !isempty(methods(f));

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

function print_info(strs...)
    global log_level
    msg = " "
    for str in strs
        msg = msg*str
    end
    log_level > 0 && println(msg)
end

function print_warn(strs...)
    print_info("\x1b[1m\x1b[31m(!)\x1b[0m ", strs...)
end

end
