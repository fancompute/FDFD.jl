export ModulatedDevice
export solve, setup_Δϵᵣ!

mutable struct ModulatedDevice <: AbstractDevice
    grid::Grid
    ϵᵣ::Array{Complex}
    Δϵᵣ::Array{Complex}
    src::Array{Complex}
    ω::AbstractVector{Float}
    Ω::Float
    nsidebands::Integer
end

function ModulatedDevice(grid::Grid, ω::AbstractVector{Float})
    ϵᵣ = ones(Complex, size(grid));
    src = zeros(Complex, size(grid));
    return Device(grid, ϵᵣ, src, ω)
end

setup_Δϵᵣ!(d::ModulatedDevice, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.Δϵᵣ, d.grid, shapes)
setup_Δϵᵣ!(d::ModulatedDevice, region, value) = _mask_values!(d.Δϵᵣ, d.grid, region, value)

function solve(d::ModulatedDevice)
    nsidebands = d.nsidebands;
    nfrequencies = 2*nsidebands+1;
    n = -nsidebands:1:nsidebands; 
    Ω = d.Ω;
    ωn = d.ω + Ω*n; 

    println("# Solver: Setup");
    println("# Solver: Number of sidebands: ", @sprintf("%d", nsidebands));
    println("# Solver: Number of frequency components: ", @sprintf("%d", nfrequencies));
    println("# Solver: Number of unknowns: ", @sprintf("%.2E", nfrequencies*length(g)));

    Ez = zeros(Complex128, nfrequencies, size(g, 1), size(g, 2));
    Hx = zeros(Complex128, nfrequencies, size(g, 1), size(g, 2));
    Hy = zeros(Complex128, nfrequencies, size(g, 1), size(g, 2));

    Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]); #TODO: check reshape vs [:]
    TΔϵ = spdiagm(ϵ₀*d.Δϵᵣ[:]); #TODO: check reshape vs [:]

    # Construct derivates
    δxb = δ(DirectionX, Backward, d.grid);
    δxf = δ(DirectionX, Forward,  d.grid);
    δyb = δ(DirectionY, Backward, d.grid);
    δyf = δ(DirectionY, Forward,  d.grid);

    # Reshape Mz into a vector
    b0 = 1im*ω*d.src[:]; #TODO: check reshape vs [:]
    b = zeros(Complex128, length(g)*nfrequencies, 1); 
    b[(nsidebands*length(g))+1:(nsidebands+1)*length(g), 1] = b0; 

    println("# Solver: Assembly");
    println("# Solver:   Calculating system matrix");
    As = Array{SparseMatrixCSC}(nfrequencies);
    Sxf = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    Sxb = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    Syf = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    Syb = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    if sharedpml
        println("# Solver:   using shared PML");
        (Sxf[1], Sxb[1], Syf[1], Syb[1]) = S_create(d.grid, ω);
        A1 = Sxb[1]*δxb*1/μ₀*Sxf[1]*δxf + Syb[1]*δyb*1/μ₀*Syf[1]*δyf;
        for i = 1:nfrequencies
            As[i] = A1 + ωn[i]^2*Tϵ;
        end
    else
        println("# Solver:   using frequency-by-frequency PML");
        for i = 1:nfrequencies
            (Sxf[i], Sxb[i], Syf[i], Syb[i]) = S_create(d.grid, ωn[i]);
            As[i] = Sxb[i]*δxb*1/μ₀*Sxf[i]*δxf + Syb[i]*δyb*1/μ₀*Syf[i]*δyf + ω[i]^2*Tϵ;
        end
    end
    if nsidebands > 0
        println("# Solver:   Calculating coupling matrix");
        stencil_Cp = spdiagm([0.5*ωn[1:end-1].^2], [1], nfrequencies, nfrequencies);
        Cp = kron(stencil_Cp, conj(TΔϵ));
        stencil_Cm = spdiagm([0.5*ωn[2:end].^2],  [-1], nfrequencies, nfrequencies);
        Cm = kron(stencil_Cm, TΔϵ);

        println("# Solver:   Calculating total matrix");
        A = blkdiag(As...) + Cp + Cm;
    else
        A = As[1]; 
    end

    println("# Solver: Solving");
    ez = dolinearsolve(A, b, matrixtype=Pardiso.COMPLEX_NONSYM);

    println("# Solver: Extracting results");
    for i = 1:nfrequencies
        Ez[i, :, :] = reshape(ez[(i-1)*length(g)+1:i*length(g)], size(g));
        Hx[i, :, :] = reshape(-1/1im/ωn[i]/μ₀*Syf[sharedpml ? 1 : i]*δyf*ez[(i-1)*length(g)+1:i*length(g)], size(g)); 
        Hy[i, :, :] = reshape(1/1im/ωn[i]/μ₀*Sxf[sharedpml ? 1 : i]*δxf*ez[(i-1)*length(g)+1:i*length(g)], size(g)); 
    end
    
    return (Ez, Hx, Hy, ωn)
end