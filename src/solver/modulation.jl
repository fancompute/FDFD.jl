export ModulatedDevice
export solve, setup_Δϵᵣ!

mutable struct ModulatedDevice{D} <: AbstractDevice{D}
    grid::Grid{D}
    ϵᵣ::Array{Complex}
    Δϵᵣ::Array{Complex}
    src::Array{Complex}
    ω::AbstractVector{Float}
    Ω::Float
    nsidebands::Integer
    sharedpml::Bool
end

"    ModulatedDevice(grid::Grid, ω::AbstractVector{Float}, Ω::Float, nsidebands::Integer; sharedpml=true)"
function ModulatedDevice(grid::Grid, ω::AbstractVector{Float}, Ω::Float, nsidebands::Integer; sharedpml=true)
    ϵᵣ = ones(Complex, size(grid));
    Δϵᵣ = zeros(Complex, size(grid));
    src = zeros(Complex, size(grid));
    return ModulatedDevice(grid, ϵᵣ, Δϵᵣ, src, ω, Ω, nsidebands, sharedpml)
end

"    ModulatedDevice(grid::Grid, ω::Float, Ω::Float, nsidebands::Integer)"
ModulatedDevice(grid::Grid, ω::Float, Ω::Float, nsidebands::Integer) = ModulatedDevice(grid, [ω], Ω, nsidebands);

"    setup_Δϵᵣ!(d::ModulatedDevice, shapes::AbstractVector{<:Shape})"
setup_Δϵᵣ!(d::ModulatedDevice, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.Δϵᵣ, d.grid, shapes)

"    setup_Δϵᵣ!(d::ModulatedDevice, region, value)"
setup_Δϵᵣ!(d::ModulatedDevice, region, value) = _mask_values!(d.Δϵᵣ, d.grid, region, value)

"    solve(d::ModulatedDevice)"
function solve(d::ModulatedDevice)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);

    ω = d.ω[1];
    nsidebands = d.nsidebands;
    nfrequencies = 2*nsidebands+1;
    n = -nsidebands:1:nsidebands; 
    Ω = d.Ω;
    ωn = ω + Ω*n;

    Tϵ = spdiagm(ϵ₀*d.ϵᵣ[:]); #TODO: check reshape vs [:]
    TΔϵ = spdiagm(ϵ₀*d.Δϵᵣ[:]); #TODO: check reshape vs [:]

    # Construct derivates
    δxb = δ(DirectionX, Backward, d.grid);
    δxf = δ(DirectionX, Forward,  d.grid);
    δyb = δ(DirectionY, Backward, d.grid);
    δyf = δ(DirectionY, Forward,  d.grid);

    # Reshape Mz into a vector
    b0 = 1im*ω*d.src[:]; #TODO: check reshape vs [:]
    b = zeros(Complex128, length(d.grid)*nfrequencies, 1); 
    b[(nsidebands*length(d.grid))+1:(nsidebands+1)*length(d.grid), 1] = b0; 

    info(FDFD.logger, "Calculating system matrix");
    info(FDFD.logger, @sprintf("Number of sidebands (frequencies): %d (%d)", nsidebands, nfrequencies));
    
    As = Array{SparseMatrixCSC}(nfrequencies);
    Sxf = Array{SparseMatrixCSC}(d.sharedpml ? 1 : nfrequencies);
    Sxb = Array{SparseMatrixCSC}(d.sharedpml ? 1 : nfrequencies);
    Syf = Array{SparseMatrixCSC}(d.sharedpml ? 1 : nfrequencies);
    Syb = Array{SparseMatrixCSC}(d.sharedpml ? 1 : nfrequencies);
    if d.sharedpml
        info(FDFD.logger, "Using shared PML");
        (Sxf[1], Sxb[1], Syf[1], Syb[1]) = S_create(d.grid, ω);
        A1 = Sxb[1]*δxb*1/μ₀*Sxf[1]*δxf + Syb[1]*δyb*1/μ₀*Syf[1]*δyf;
        for i = 1:nfrequencies
            As[i] = A1 + ωn[i]^2*Tϵ;
        end
    else
        info(FDFD.logger, "Using frequency-by-frequency PML");
        for i = 1:nfrequencies
            (Sxf[i], Sxb[i], Syf[i], Syb[i]) = S_create(d.grid, ωn[i]);
            As[i] = Sxb[i]*δxb*1/μ₀*Sxf[i]*δxf + Syb[i]*δyb*1/μ₀*Syf[i]*δyf + ωn[i]^2*Tϵ;
        end
    end
    if nsidebands > 0
        info(FDFD.logger, "Calculating coupling matrix");
        stencil_Cp = spdiagm([0.5*ωn[1:end-1].^2], [1], nfrequencies, nfrequencies);
        Cp = kron(stencil_Cp, conj(TΔϵ));
        stencil_Cm = spdiagm([0.5*ωn[2:end].^2],  [-1], nfrequencies, nfrequencies);
        Cm = kron(stencil_Cm, TΔϵ);

        info(FDFD.logger, "Calculating total matrix");
        A = blkdiag(As...) + Cp + Cm;
    else
        A = As[1]; 
    end

    ez = dolinearsolve(A, b, matrixtype=Pardiso.COMPLEX_NONSYM);

    info(FDFD.logger, "Extracting results");

    fields = Array{FieldTM}(nfrequencies);
    for i = 1:nfrequencies
        ezi = ez[(i-1)*length(d.grid)+1:i*length(d.grid)];
        hxi = -1/1im/ωn[i]/μ₀*Syf[d.sharedpml ? 1 : i]*δyf*ezi; 
        hyi = 1/1im/ωn[i]/μ₀*Sxf[d.sharedpml ? 1 : i]*δxf*ezi; 
        fields[i] = FieldTM(d.grid, ezi, hxi, hyi)
    end
    
    return (fields, ωn)
end