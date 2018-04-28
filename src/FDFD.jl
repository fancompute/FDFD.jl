module FDFD
using Memento
const logger = getlogger(current_module());
__init__() = Memento.register(logger);

using GeometryPrimitives, StaticArrays, Pardiso

include("./types.jl");
include("./grid.jl");
include("./data.jl");
include("./device.jl");
include("./pml.jl");
include("./solver/driven.jl");
include("./solver/eigen.jl");
include("./solver/modulation.jl");
include("./solver/nonlinear.jl");
include("./plot.jl");

export poynting, flux_surface, unwrap

"""
    δ(w::Direction, s::DerivativeDirection, g::Grid)

Compute the derivative operators for a given geometry. w should be one of "x" or "y"
indicating the direction of the derivative. s should be one of "f" or "b" to indicate
forward or backward derivative, respectively.
"""
function δ(w::Direction, s::DerivativeDirection, g::Grid{K}) where {K}
    if K == 2
	    (Nx, Ny) = size(g);
    else 
        Nx = size(g, 1);
        Ny = 1;
    end

    if w == DirectionX
        if s == Forward
            δxf = 1/dx(g)*spdiagm([ones(Nx-1), -ones(Nx), 1], [1, 0, -Nx+1]);
            return kron(speye(Ny), δxf)
        else
            δxb = 1/dx(g)*spdiagm([-ones(Nx-1), ones(Nx), -1], [-1, 0, Nx-1]);
            return kron(speye(Ny), δxb)
        end
    end
    if w == DirectionY
        if s == Backward
            δyf = 1/dy(g)*spdiagm([ones(Ny-1), -ones(Ny), 1], [1, 0, -Ny+1]);
            return kron(δyf, speye(Nx))
        else
            δyb = 1/dy(g)*spdiagm([-ones(Ny-1), ones(Ny), -1], [-1, 0, Ny-1]);
            return kron(δyb, speye(Nx))
        end
    end
end


"""
    grid_average(center_array, w)


"""
function grid_average(centerarray::AbstractArray, w::Direction)
    ndims(centerarray) == 1 && return (centerarray+circshift(centerarray, (1)))/2
    w == DirectionX && return (centerarray+circshift(centerarray, (1, 0)))/2;
    w == DirectionY && return (centerarray+circshift(centerarray, (0, 1)))/2;
    return centerarray
end


"""
    poynting(Ez, Hx, Hy)


"""
function poynting(field::Field)
    if isa(field, FieldTM)
        Ez_x = grid_average(field.Ez, DirectionX);
        Ez_y = grid_average(field.Ez, DirectionY);
        Sx = -0.5*real.(Ez_x.*conj(field.Hy));
        Sy = 0.5*real.(Ez_y.*conj(field.Hx));
        return Flux2D(field.grid, Sx, Sy)
    end
    if isa(field, FieldTE)
        error("Not implemented yet...");
    end
    error("Invalid polarization");
end

"   flux_surface(field::Field, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction) = flux_surface(poynting(field), ptmid, width, nrm)"
flux_surface(field::Field, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction) = flux_surface(poynting(field), ptmid, width, nrm);

"    flux_surface(flux::Flux2D, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction)"
function flux_surface(flux::Flux2D, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction)
    (x0, y0) = coord2ind(flux.grid, ptmid);
    if nrm == DirectionX
        indx = x0;
        if isinf(width)
            indy = 1:flux.grid.N[2];
        else
            y1 = y2ind(flux.grid,ptmid[2]-width/2);
            y2 = y2ind(flux.grid,ptmid[2]+width/2);
            indy = y1:1:y2;
        end
        return sum(flux.Sx[indx,indy])*dy(flux.grid)
    elseif nrm == DirectionY
        error("Not implemented yet...")
    end
end

function dolinearsolve(A::SparseMatrixCSC, b::Array; matrixtype=Pardiso.COMPLEX_NONSYM)
    tic();
    info(FDFD.logger, "Performing linear solve");
    info(FDFD.logger, @sprintf("Problem unknowns: %.2E", length(b)));
    pardiso_success = false;
    try
        ps = PardisoSolver();
        set_matrixtype!(ps, matrixtype);
        set_solver!(ps, Pardiso.DIRECT_SOLVER);
        pardisoinit(ps);
        x = solve(ps, A, b);
        pardiso_success = true;
    catch
        debug(FDFD.logger, "Pardiso solver has failed, falling back to lufact()");
    end
    if ~pardiso_success
        x = lufact(A)\b;
    end
    info(FDFD.logger, @sprintf("Time to solve: %.2f minutes", toq()/60));
    return x
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

end