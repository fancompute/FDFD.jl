module fdfd

const epsilon0 = 8.85418782e-12;
const mu0 = 1.25663706e-6;
const c0 = sqrt(1/epsilon0/mu0);
const eta0 = sqrt(mu0/epsilon0);

include("./helpers.jl");
include("./datastructs.jl");
include("./pml.jl");
include("./eigen.jl");
include("./driven.jl");
include("./modulation.jl");
include("./plot.jl");

export assign_src!, assign_src_func!, assign_src_point!, assign_src_mode!, assign_epsr!, assign_epsr_func!
export coord2ind, poyntingTM

function dws(w, s, geom)
    if typeof(geom) == Geometry
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
            dxf = 1/dx(geom)*spdiagm([ones(Nx-1), -ones(Nx), 1], [1, 0, -Nx+1]);
            return kron(speye(Ny), dxf)
        else
            dxb = 1/dx(geom)*spdiagm([-ones(Nx-1), ones(Nx), -1], [-1, 0, Nx-1]);
            return kron(speye(Ny), dxb)
        end
    else
        if s == "f"
            dyf = 1/dy(geom)*spdiagm([ones(Ny-1), -ones(Ny), 1], [1, 0, -Ny+1]);
            return kron(dyf, speye(Nx))
        else
            dyb = 1/dy(geom)*spdiagm([-ones(Ny-1), ones(Ny), -1], [-1, 0, Ny-1]);
            return kron(dyb, speye(Nx))
        end
    end
end

function grid_average(center_array, w)
    ndims(center_array) == 1 && return (center_array+circshift(center_array, (1)))/2
    w == "x" && return (center_array+circshift(center_array, (1, 0)))/2;
    w == "y" && return (center_array+circshift(center_array, (0, 1)))/2;
    return center_array
end

function assign_src!(geom::Geometry, region, value)
    mask = [region(x, y) for x in xc(geom), y in yc(geom)];
    geom.src[mask] = value;
end

function assign_src_point!(geom::Geometry, xy)
    (indx, indy) = coord2ind(geom, xy);
    geom.src[indx, indy] = 1im;
end

function assign_src_func!(geom::Geometry, region, value_func)
    mask = [region(x, y) for x in xc(geom), y in yc(geom)];
    value_computed = [value_func(x, y) for x in xc(geom), y in yc(geom)];
    geom.src[mask] = value_computed[mask];
end

function assign_src_mode!(geom::Geometry, pol, omega, beta_est, src_xy, src_normal, Nsrc)
    (src_ind_x, src_ind_y) = coord2ind(geom, src_xy);

    NN = Int64(round((Nsrc-1)/2))
    Nsrc = 2*NN+1
    if src_normal == "x"
        inds_x = src_ind_x;
        inds_y = src_ind_y+(-NN:NN);
        dh = dy(geom);
    elseif src_normal == "y"
        inds_x = src_ind_x;
        inds_y = src_ind_y+(-NN:NN);
        dh = dx(geom);
    else
        error("Invalid src_normal value. Must be x or y")
    end

    epsr_src = geom.epsr[inds_x, inds_y];
    src_range = (0.0, Nsrc*dh);
    geom1D = Geometry1D(src_range, epsr_src);
    (beta, output_vector) = solve_eigen_1D(geom1D, pol, omega, beta_est, 1);
    geom.src[inds_x, inds_y] = 1im*output_vector;
end

function assign_epsr!(geom::Geometry, region, value)
    mask = [region(x, y) for x in xc(geom), y in yc(geom)];
    geom.epsr[mask] = value;
end

function assign_epsr_func!(geom::Geometry, region, value_func)
    mask = [region(x, y) for x in xc(geom), y in yc(geom)];
    value_computed = [value_func(x, y) for x in xc(geom), y in yc(geom)];
    geom.epsr[mask] = value_computed[mask];
end

function poyntingTM(Ez, Hx, Hy)
    Ez_x = grid_average(Ez, "x")
    Sx = -1/2*real(Ez_x.*conj(Hy));
    Ez_y = grid_average(Ez, "y");
    Sy = 1/2*real(Ez_y.*conj(Hx));
    return (Sx, Sy)
end

end