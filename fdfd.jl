module fdfd

const ϵ₀ = 8.85418782e-12;
const μ₀ = 1.25663706e-6;
const c₀ = sqrt(1/ϵ₀/μ₀);
const η₀ = sqrt(μ₀/ϵ₀);

include("./helpers.jl");
include("./datastructs.jl");
include("./pml.jl");
include("./eigen.jl");
include("./driven.jl");
include("./modulation.jl");
include("./plot.jl");

export assign_src!, assign_src_func!, assign_src_point!, assign_src_mode!, assign_epsr!, assign_epsr_func!
export coord2ind, poyntingTM, flux_direction
export ϵ₀, μ₀, c₀, η₀

"""
    δ(w, s, geom::Geometry)

Compute the derivative operators for a given geometry. w should be one of "x" or "y"
indicating the direction of the derivative. s should be one of "f" or "b" to indicate
forward or backward derivative, respectively.
"""
function δ(w, s, geom::Geometry)
    if typeof(geom) == Geometry2D
	    Nx = geom.N[1];
	    Ny = geom.N[2];
    elseif typeof(geom) == Geometry1D
        Nx = geom.N;
        Ny = 1;
    else
        error("Unkown geometry type!")
    end

    if w == "x"
        if s == "f"
            δxf = 1/dx(geom)*spdiagm([ones(Nx-1), -ones(Nx), 1], [1, 0, -Nx+1]);
            return kron(speye(Ny), δxf)
        else
            δxb = 1/dx(geom)*spdiagm([-ones(Nx-1), ones(Nx), -1], [-1, 0, Nx-1]);
            return kron(speye(Ny), δxb)
        end
    else
        if s == "f"
            δyf = 1/dy(geom)*spdiagm([ones(Ny-1), -ones(Ny), 1], [1, 0, -Ny+1]);
            return kron(δyf, speye(Nx))
        else
            δyb = 1/dy(geom)*spdiagm([-ones(Ny-1), ones(Ny), -1], [-1, 0, Ny-1]);
            return kron(δyb, speye(Nx))
        end
    end
end


"""
    grid_average(center_array, w)


"""
function grid_average(centerarray, w)
    ndims(centerarray) == 1 && return (centerarray+circshift(centerarray, (1)))/2
    w == "x" && return (centerarray+circshift(centerarray, (1, 0)))/2;
    w == "y" && return (centerarray+circshift(centerarray, (0, 1)))/2;
    return centerarray
end


"""
    assign_src!(geom::Geometry2D, region, value)


"""
function assign_src!(geom::Geometry2D, region, value)
    mask = [region(x, y) for x in xc(geom), y in yc(geom)];
    if iscallable(value)
        value_assigned = [value(x, y) for x in xc(geom), y in yc(geom)];
    else
        value_assigned = value;
    end
    geom.src[mask] = value_assigned;
end


"""
    assign_src_point!(geom::Geometry2D, xy)


"""
function assign_src_point!(geom::Geometry2D, xy)
    (indx, indy) = coord2ind(geom, xy);
    geom.src[indx, indy] = 1im;
end


"""
    assign_src_mode!(geom::Geometry2D, pol, omega, beta_est, src_xy, src_normal, Nsrc)


"""
function assign_src_mode!(geom::Geometry2D, polarization, ω, estimatedβ, srcxy, srcnormal, srcpoints)
    (indx, indy) = coord2ind(geom, srcxy);

    M = Int64(round((srcpoints-1)/2))
    srcpoints = 2*M+1
    if src_normal == "x"
        indx = indx;
        indy = indy+(-M:M);
        dh = dy(geom);
    elseif src_normal == "y"
        indx = indx;
        indy = indy+(-M:M);
        dh = dx(geom);
    else
        error("Invalid src_normal value. Must be x or y")
    end

    ϵᵣ = geom.ϵᵣ[indx, indy];
    srange = (0.0, Nsrc*dh);
    geom1D = Geometry1D(srange, ϵᵣ);
    (β, vector) = solve_eigen_1D(geom1D, polarization, ω, estimatedβ, 1);
    geom.src[indx, indy] = 1im*vector;
end


"""
    assign_epsr!(geom::Geometry, region, value)


"""
function assign_epsr!(geom::Geometry, region, value)
    if typeof(geom) == Geometry2D
        mask = [region(x, y) for x in xc(geom), y in yc(geom)];
        if iscallable(value)
            value_assigned = [value(x, y) for x in xc(geom), y in yc(geom)];
        else
            value_assigned = value;
        end
    elseif typeof(geom) == Geometry1D
        mask = [region(x) for x in xc(geom)];
        if iscallable(value)
            value_assigned = [value(x) for x in xc(geom)];
        else
            value_assigned = value;
        end
    else
        error("Unkown geometry type!")
    end
    geom.ϵᵣ[mask] = value_assigned;
end


"""
    poyntingTM(Ez, Hx, Hy)


"""
function poyntingTM(Ez, Hx, Hy)
    Ez_x = grid_average(Ez, "x")
    Sx = -1/2*real(Ez_x.*conj(Hy));
    Ez_y = grid_average(Ez, "y");
    Sy = 1/2*real(Ez_y.*conj(Hx));
    return (Sx, Sy)
end

function flux_direction(dir_normal, pt1, pt2, geom, Ez, Hx, Hy)
    if dir_normal == "x"
        (ind0, _) = coord2ind(geom, (pt1, geom.yrange[1]));
        (ind1, _) = coord2ind(geom, (pt2, geom.yrange[1]));
        dh = dy(geom);
    elseif dir_normal == "y"
        error("Not implemented yet...")
    end
    ind_cells = ind0:ind1;
    N_cells = ind1-ind0+1;

    N_freqs = size(Ez)[1];
    Px = zeros(Real,N_freqs,N_cells);

    for i = 1:N_freqs
        (Sx, _) = poyntingTM(Ez[i,:,:], Hx[i,:,:], Hy[i,:,:]);
        Px[i, :] = sum(Sx[ind0:ind1,:],2)*dh;
    end

    x_coords = xc(geom);
    x_coords = x_coords[ind0:ind1];
    Px = Px.';
    return (x_coords, Px)
end

end