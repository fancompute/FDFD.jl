module fdfd

const epsilon0 = 8.85418782e-12;
const mu0 = 1.25663706e-6;
const c0 = sqrt(1/epsilon0/mu0);
const eta0 = sqrt(mu0/epsilon0);

include("./pml.jl");
include("./eigen.jl");
include("./driven.jl");
include("./modulation.jl");

export assign_val!, assign_val_func!, assign_modal_source!
export coord2ind, poyntingTM

function dws(w, s, N, xrange, yrange)
	Nx = N[1];
	Ny = N[2];
    dx = (xrange[2]-xrange[1])/Nx;
    dy = (yrange[2]-yrange[1])/Ny;

    if w == "x"
        if s == "f"
            dxf = 1/dx*spdiagm([ones(Nx-1), -ones(Nx), 1], [1, 0, -Nx+1]);
            return kron(speye(Ny), dxf)
        else
            dxb = 1/dx*spdiagm([-ones(Nx-1), ones(Nx), -1], [-1, 0, Nx-1]);
            return kron(speye(Ny), dxb)
        end
    else
        if s == "f"
            dyf = 1/dy*spdiagm([ones(Ny-1), -ones(Ny), 1], [1, 0, -Ny+1]);
            return kron(dyf, speye(Nx))
        else
            dyb = 1/dy*spdiagm([-ones(Ny-1), ones(Ny), -1], [-1, 0, Ny-1]);
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

function assign_val!(val_array, region, value, xrange, yrange)
    (Nx, Ny) = size(val_array);
	dx = (xrange[2]-xrange[1])/Nx;
    dy = (yrange[2]-yrange[1])/Ny;

    xc = xrange[1]+dx*(0.5:1:Nx);
    yc = yrange[1]+dy*(0.5:1:Ny);

    mask = [region(x, y) for x in xc, y in yc];
    val_array[mask] = value;
end

function assign_val_func!(val_array, region, value_func, xrange, yrange)
    (Nx, Ny) = size(val_array);
    dx = (xrange[2]-xrange[1])/Nx;
    dy = (yrange[2]-yrange[1])/Ny;

    xc = xrange[1]+dx*(0.5:1:Nx);
    yc = yrange[1]+dy*(0.5:1:Ny);

    mask = [region(x, y) for x in xc, y in yc];
    value_computed = [value_func(x, y) for x in xc, y in yc];
    val_array[mask] = value_computed[mask];
end

function assign_modal_source!(pol, omega, beta_est, Jz, eps_r, src_xy, src_normal, Nsrc, N, xrange, yrange)
    (src_ind_x, src_ind_y) = coord2ind(src_xy, xrange, yrange, N);
    dx = (xrange[2]-xrange[1])/N[1];
    dy = (yrange[2]-yrange[1])/N[2];

    NN = Int64(round((Nsrc-1)/2))
    Nsrc = 2*NN+1
    if src_normal == "x"
        inds_x = src_ind_x;
        inds_y = src_ind_y+(-NN:NN);
        dh = dy;
    elseif src_normal == "y"
        inds_x = src_ind_x;
        inds_y = src_ind_y+(-NN:NN);
        dh = dx;
    else
        error("Invalid src_normal value. Must be x or y")
    end

    eps_r_src = eps_r[inds_x, inds_y];
    src_range = (0, Nsrc*dh);
    (beta, output_vector) = solve_eigen_1D(pol, omega, beta_est, 1, src_range, eps_r_src);
    Jz[inds_x, inds_y] = 1im*output_vector;
end

function coord2ind(xy, xrange, yrange, N)
    indx = Int64(round((xy[1]-xrange[1])/(xrange[2]-xrange[1])*N[1])+1); 
    indy = Int64(round((xy[2]-yrange[1])/(yrange[2]-yrange[1])*N[2])+1);
    return (indx, indy)
end

function poyntingTM(Ez, Hx, Hy)
    Ez_x = grid_average(Ez, "x")
    Sx = -1/2*real(Ez_x.*conj(Hy));
    Ez_y = grid_average(Ez, "y");
    Sy = 1/2*real(Ez_y.*conj(Hx));
    return (Sx, Sy)
end

end