module FDFD

using GeometryPrimitives, StaticArrays, Pardiso, Memento

const logger = getlogger(current_module());
__init__() = Memento.register(logger);

const ϵ₀ = 8.85418782e-12;
const μ₀ = 1.25663706e-6;
const c₀ = sqrt(1/ϵ₀/μ₀);
const η₀ = sqrt(μ₀/ϵ₀);

include("./types.jl");
include("./grid.jl");
include("./device.jl");
include("./pml.jl");
include("./solver/driven.jl");
include("./solver/eigen.jl");
include("./solver/modulation.jl");
include("./solver/nonlinear.jl");
#include("./plot.jl");

export ϵ₀, μ₀, c₀, η₀
export poynting, flux_direction, unwrap

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
function poynting(polarization::Polarization, Ez, Hx, Hy)
    if polarization == TM
        Ez_x = grid_average(Ez, DirectionX);
        Ez_y = grid_average(Ez, DirectionY);
        Sx = -1/2*real(Ez_x.*conj(Hy));
        Sy = 1/2*real(Ez_y.*conj(Hx));
    elseif polarization == TE
        error("Not implemented yet...");
    else
        error("Invalid polarization");
    end
    return (Sx, Sy)
end

function flux_direction(dir_normal::Direction, pt1, pt2, geom, Ez, Hx, Hy)
    if dir_normal == DirectionX
        (ind0, _) = coord2ind(geom, (pt1, geom.yrange[1]));
        (ind1, _) = coord2ind(geom, (pt2, geom.yrange[1]));
        dh = dy(geom);
    elseif dir_normal == DirectionY
        error("Not implemented yet...")
    end
    ind_cells = ind0:ind1;
    N_cells = ind1-ind0+1;

    N_freqs = size(Ez)[1];
    Px = zeros(Real,N_freqs,N_cells);

    for i = 1:N_freqs
        (Sx, _) = poynting(TM, Ez[i,:,:], Hx[i,:,:], Hy[i,:,:]);
        Px[i, :] = sum(Sx[ind0:ind1,:],2)*dh;
    end

    x_coords = xc(geom);
    x_coords = x_coords[ind0:ind1];
    Px = Px.';
    return (x_coords, Px)
end

function dolinearsolve(A, b; matrixtype=Pardiso.COMPLEX_NONSYM)
    global warned_pardiso
    tic();
    pardiso_success = false;
    try
        ps = PardisoSolver();
        set_matrixtype!(ps, matrixtype);
        set_solver!(ps, Pardiso.DIRECT_SOLVER);
        pardisoinit(ps);
        x = solve(ps, A, b);
        pardiso_success = true;
    catch
        #warn(logger, "Pardiso solver has failed, falling back to lufact()");
    end
    if ~pardiso_success
        x = lufact(A)\b;
    end
    info(logger, @sprintf("Solve completed in %.2f min", toq()/60));
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