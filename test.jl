include("./fdfd.jl");
using fdfd

omega = 2*pi*200e12
N = (10, 10);
Npml = (0,0)
xrange = (-1e-6, 1e-6);
yrange = (-1e-6, 1e-6);

eps_r = ones(N);
region_func(x,y) = -0.25e-6<=x<=0.25e-6 && -0.1e-6<=y<=0.1e-6;
fdfd.assign_val!(eps_r, region_func, 2, xrange, yrange);

Jz = zeros(N)
Jz[5,5]=1

(A,b)=fdfd.solve_TM(omega, xrange, yrange, eps_r, Jz, Npml)
spy(A)
