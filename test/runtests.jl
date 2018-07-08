using FDFD, GeometryPrimitives
using Base.Test

set_log_level!(0);
ENV["PARDISOLICMESSAGE"]=0

function compare_solvers_dipole()
    grd = Grid(0.01, [15 15], [-3 3], [-3 3]);
    dev = Device(grd, 2*π*200e12);
    setup_src!(dev, Point(0, 0));

    ENV["FDFD_SOLVER"] = "julia"
    field_julia = solve(dev);
    ENV["FDFD_SOLVER"] = "pardiso"
    field_pardiso = solve(dev);
    ENV["FDFD_SOLVER"] = "julia"

    equality_real = real.(field_julia) ≈ real.(field_pardiso)
    equality_imag = imag.(field_julia) ≈ imag.(field_pardiso)
    return equality_real && equality_imag
end

function compare_solvers_wg()
    grd = Grid(0.02, [15, 15], [0.0, 10.0], [-1.0, 1.0]);
    dev = Device(grd, 2*π*200e12);
    setup_ϵᵣ!(dev, [Box([5.0, 0.0, 0.0], [Inf, 0.3, Inf], eye(3), 12)])
    add_mode!(dev, Mode(TM, x̂, 3.5, Point(1.0, 0), 0.8));

    ENV["FDFD_SOLVER"] = "julia"
    field_julia = solve(dev);
    ENV["FDFD_SOLVER"] = "pardiso"
    field_pardiso = solve(dev);
    ENV["FDFD_SOLVER"] = "julia"

    equality_real = real.(field_julia) ≈ real.(field_pardiso)
    equality_imag = imag.(field_julia) ≈ imag.(field_pardiso)
    return equality_real && equality_imag
end

function compare_solvers_modwg()
    nsidebands = 1
    L = 5.0
    grd = Grid(0.01, [15 15], [0.0 L], [-1.0 1.0]);
    dev = ModulatedDevice(grd, 2π*1.939e14, 4.541e14, nsidebands)
    a = 0.2202
    q = 2.9263
    region_wg(x,y) = -a/2<=y<=a/2
    region_modulation(x,y) = 1<=x<=(L-1) && -a/2<=y<=0
    modulation_function(x,y) = exp.(1im*q*x)
    setup_ϵᵣ!(dev, region_wg, 12.25)
    setup_Δϵᵣ!(dev, region_modulation, modulation_function)
    add_mode!(dev, Mode(TM, x̂, 3.5, [0.2, 0], 4*a))

    ENV["FDFD_SOLVER"] = "julia"
    field_julia = solve(dev);
    ENV["FDFD_SOLVER"] = "pardiso"
    field_pardiso = solve(dev);
    ENV["FDFD_SOLVER"] = "julia"

    equality_real = all([real.(field_julia[i]) ≈ real.(field_pardiso[i]) for i = 1:2*nsidebands+1])
    equality_imag = all([imag.(field_julia[i]) ≈ imag.(field_pardiso[i]) for i = 1:2*nsidebands+1])
    return equality_real && equality_imag
end

@testset "FDFD" begin
    @testset "Solver" begin
        @test compare_solvers_dipole()
        @test compare_solvers_wg()
        @test compare_solvers_modwg()
    end
end
