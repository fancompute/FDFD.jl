module fdfd

const epsilon0 = 8.85418782e-12;
const mu0 = 1.25663706e-6;
const c0 = sqrt(1/epsilon0/mu0);
const eta0 = sqrt(mu0/epsilon0);

include("./pml.jl");

export assign_val, solve_TM

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

    val_array[[region(x,y) for x in xc, y in yc]] = value;
end

function solve_TM(omega, xrange, yrange, eps_r, Jz, Npml)
    N = size(eps_r);

    T_eps_z = spdiagm(epsilon0*eps_r[:]);

    jz = Jz[:];

    Hx = zeros(Complex128, N);
    Hy = zeros(Complex128, N);
    Ez = zeros(Complex128, N);

    (Sxf, Sxb, Syf, Syb) = S_create(omega, N, Npml, xrange, yrange);

    # Construct derivates
    Dyb = Syb*dws("y", "b", N, xrange, yrange);
    Dxb = Sxb*dws("x", "b", N, xrange, yrange);
    Dxf = Sxf*dws("x", "f", N, xrange, yrange);
    Dyf = Syf*dws("y", "f", N, xrange, yrange);

    # Construct system matrix
    A = Dxf*mu0^-1*Dxb + Dyf*mu0^-1*Dyb + omega^2*T_eps_z;
    b = 1im*omega*jz;

    ez = A\b;

    hx = -1/1im/omega/mu0*Dyb*ez;
    hy = 1/1im/omega/mu0*Dxb*ez;

    Hx = reshape(hx, N);
    Hy = reshape(hy, N);
    Ez = reshape(ez, N);

    return (Ez, Hx, Hy)
end

end