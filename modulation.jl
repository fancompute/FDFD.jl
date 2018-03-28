export Modulator
export solve_modulation_TM
export setΔϵᵣ!

using Pardiso

mutable struct Modulator
    geom::Geometry2D
    Ω₁::Real
    Ω₂::Real
    nsidebands::Integer
    Δϵᵣ::Array{Complex,2}

    function Modulator(geom::Geometry2D, Ω::Real, nsidebands::Integer)
        return Modulator(geom, Ω, 0, nsidebands)
    end
    function Modulator(geom::Geometry2D, Ω₁::Real, Ω₂::Real, nsidebands::Integer)
        Δϵᵣ = zeros(Complex128, geom.N)
        return new(geom, Ω₁, Ω₂, nsidebands, Δϵᵣ)
    end
end

function setΔϵᵣ!(mod::Modulator, region, value)
    mask = [region(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
    if iscallable(value)
        value_unmasked = [value(x, y) for x in xc(mod.geom), y in yc(mod.geom)];
        value_assigned = value_unmasked[mask];
    else
        value_assigned = value;
    end
    mod.Δϵᵣ[mask] = value_assigned;
end

function solve_modulation_TM(mod::Modulator, ω₀::Real; verbose=true, sharedpml=true)
    geom = mod.geom;
    nsidebands = mod.nsidebands;
    nfrequencies = 2*nsidebands+1;
    n = -nsidebands:1:nsidebands; 
    N = geom.N;
    M = prod(N);
    Ω = mod.Ω₁;
    ω = ω₀ + Ω*n; 

    println("# Solver: Setup");
    println("# Solver: Number of sidebands: ", @sprintf("%d", nsidebands));
    println("# Solver: Number of frequency components: ", @sprintf("%d", nfrequencies));
    println("# Solver: Number of unknowns: ", @sprintf("%.2E", nfrequencies*M));

    Ez = zeros(Complex128, nfrequencies, N[1], N[2]);
    Hx = zeros(Complex128, nfrequencies, N[1], N[2]);
    Hy = zeros(Complex128, nfrequencies, N[1], N[2]);

    Tϵ = spdiagm(ϵ₀*geom.ϵᵣ[:]); #TODO: check reshape vs [:]
    TΔϵ = spdiagm(ϵ₀*mod.Δϵᵣ[:]); #TODO: check reshape vs [:]

    # Construct derivates
    δxb = δ("x", "b", geom);
    δxf = δ("x", "f", geom);
    δyb = δ("y", "b", geom);
    δyf = δ("y", "f", geom); 

    # Reshape Mz into a vector
    b0 = 1im*ω₀*geom.src[:]; #TODO: check reshape vs [:]
    b = zeros(Complex128, M*nfrequencies,1); 
    b[(nsidebands*M)+1:(nsidebands+1)*M,1] = b0; 

    println("# Solver: Assembly");
    println("# Solver:   Calculating system matrix");
    As = Array{SparseMatrixCSC}(nfrequencies);
    Sxf = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    Sxb = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    Syf = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    Syb = Array{SparseMatrixCSC}(sharedpml ? 1 : nfrequencies);
    if sharedpml
        println("# Solver:   using shared PML");
        (Sxf[1], Sxb[1], Syf[1], Syb[1]) = S_create(ω₀, N, geom.Npml, geom.xrange, geom.yrange);
        A1 = Sxb[1]*δxb*1/μ₀*Sxf[1]*δxf + Syb[1]*δyb*1/μ₀*Syf[1]*δyf;
        for i = 1:nfrequencies
            As[i] = A1 + ω[i]^2*Tϵ;
        end
    else
        println("# Solver:   using frequency-by-frequency PML");
        for i = 1:nfrequencies
            (Sxf[i], Sxb[i], Syf[i], Syb[i]) = S_create(ω[i], N, geom.Npml, geom.xrange, geom.yrange);
            As[i] = Sxb[i]*δxb*1/μ₀*Sxf[i]*δxf + Syb[i]*δyb*1/μ₀*Syf[i]*δyf + ω[i]^2*Tϵ;
        end
    end
    if nsidebands > 0
        println("# Solver:   Calculating coupling matrix");
        stencil_Cp = spdiagm([0.5*ω[1:end-1].^2], [1], nfrequencies, nfrequencies);
        Cp = kron(stencil_Cp, conj(TΔϵ));
        stencil_Cm = spdiagm([0.5*ω[2:end].^2],  [-1], nfrequencies, nfrequencies);
        Cm = kron(stencil_Cm, TΔϵ);

        println("# Solver:   Calculating total matrix");
        A = blkdiag(As...) + Cp + Cm;
    else
        A = As[1]; 
    end

    println("# Solver: Solving");
    ez = dolinearsolve(A, b, matrixtype=Pardiso.COMPLEX_NONSYM, verbose=verbose);

    println("# Solver: Extracting results");
    for i = 1:nfrequencies
        Ez[i, :, :] = reshape(ez[(i-1)*M+1:i*M], N);
        Hx[i, :, :] = reshape(-1/1im/ω[i]/μ₀*Syf[sharedpml ? 1 : i]*δyf*ez[(i-1)*M+1:i*M], N); 
        Hy[i, :, :] = reshape(1/1im/ω[i]/μ₀*Sxf[sharedpml ? 1 : i]*δxf*ez[(i-1)*M+1:i*M], N); 
    end
    
    return (Ez, Hx, Hy, ω)
end