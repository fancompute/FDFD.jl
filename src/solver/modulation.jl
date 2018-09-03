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
    modes::Array{Mode}
end

"    ModulatedDevice(grid::Grid, ω::AbstractVector{Float}, Ω::Float, nsidebands::Integer; sharedpml::Bool=true)"
function ModulatedDevice(grid::Grid, ω::AbstractVector{Float}, Ω::Float, nsidebands::Integer; sharedpml::Bool=true)
    ϵᵣ = ones(Complex, size(grid));
    Δϵᵣ = zeros(Complex, size(grid));
    src = zeros(Complex, size(grid));
    return ModulatedDevice(grid, ϵᵣ, Δϵᵣ, src, ω, Ω, nsidebands, sharedpml, Array{Mode}(undef, 0))
end

"    ModulatedDevice(grid::Grid, ω::Float, Ω::Float, nsidebands::Integer; sharedpml::Bool=true)"
ModulatedDevice(grid::Grid, ω::Float, Ω::Float, nsidebands::Integer; sharedpml::Bool=true) =
 ModulatedDevice(grid, [ω], Ω, nsidebands, sharedpml=sharedpml);

"    setup_Δϵᵣ!(d::ModulatedDevice, shapes::AbstractVector{<:Shape})"
setup_Δϵᵣ!(d::ModulatedDevice, shapes::AbstractVector{<:Shape}) = _compose_shapes!(d.Δϵᵣ, d.grid, shapes)

"    setup_Δϵᵣ!(d::ModulatedDevice, region, value)"
setup_Δϵᵣ!(d::ModulatedDevice, region, value) = _mask_values!(d.Δϵᵣ, d.grid, region, value)

"    solve(d::ModulatedDevice)"
function solve(d::ModulatedDevice)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d);

    Nω = length(d.ω);
    nsidebands = d.nsidebands;
    nfrequencies = 2*nsidebands+1;
    n = -nsidebands:1:nsidebands;
    Ω = d.Ω;

    fields = Array{FieldTM}(undef, Nω, nfrequencies);

    for i in eachindex(d.ω)
        @info "Frequency: $i/$Nω"
        ω = d.ω[i];
        ωn = ω .+ Ω*n;

        length(d.modes) > 0 && ( d.src=zeros(Complex, size(d.grid)) ) # Only reset if we are using modes
        for mode in d.modes
            # TODO: handle the polarization here
            setup_mode!(d, TM, ω, mode.neff, mode.pt, mode.dir, mode.width);
        end

        Tϵ = sparse(Diagonal(ϵ₀*d.ϵᵣ[:])); #TODO: check reshape vs [:]
        TΔϵ = sparse(Diagonal(ϵ₀*d.Δϵᵣ[:])); #TODO: check reshape vs [:]

        # Construct derivates
        δxb = δ(x̂, Backward, d.grid);
        δxf = δ(x̂, Forward,  d.grid);
        δyb = δ(ŷ, Backward, d.grid);
        δyf = δ(ŷ, Forward,  d.grid);

        # Reshape Mz into a vector
        b0 = 1im*ω*d.src[:]; #TODO: check reshape vs [:]
        b = zeros(ComplexF64, length(d.grid)*nfrequencies, 1);
        b[(nsidebands*length(d.grid))+1:(nsidebands+1)*length(d.grid), 1] = b0;

        @info "Calculating: system matrix"
        @info "Sidebands (freqs): $nsidebands ($nfrequencies)"

        As = Array{SparseMatrixCSC}(undef, nfrequencies);
        Sxf = Array{SparseMatrixCSC}(undef, d.sharedpml ? 1 : nfrequencies);
        Sxb = Array{SparseMatrixCSC}(undef, d.sharedpml ? 1 : nfrequencies);
        Syf = Array{SparseMatrixCSC}(undef, d.sharedpml ? 1 : nfrequencies);
        Syb = Array{SparseMatrixCSC}(undef, d.sharedpml ? 1 : nfrequencies);
        if d.sharedpml
            @info "PML: shared"
            (Sxf[1], Sxb[1], Syf[1], Syb[1]) = S_create(d.grid, ω);
            A1 = Sxb[1]*δxb*1/μ₀*Sxf[1]*δxf + Syb[1]*δyb*1/μ₀*Syf[1]*δyf;
            for j = 1:nfrequencies
                As[j] = A1 + ωn[j]^2*Tϵ;
            end
        else
            @info "PML: per-frequency"
            for j = 1:nfrequencies
                (Sxf[j], Sxb[j], Syf[j], Syb[j]) = S_create(d.grid, ωn[j]);
                As[j] = Sxb[j]*δxb*1/μ₀*Sxf[j]*δxf + Syb[j]*δyb*1/μ₀*Syf[j]*δyf + ωn[j]^2*Tϵ;
            end
        end
        if nsidebands > 0
            @info "Calculating: coupling matrix"
            stencil_Cp = spdiagm( 1 => 0.5*ωn[1:end-1].^2);
            Cp = kron(stencil_Cp, conj(TΔϵ));
            stencil_Cm = spdiagm(-1 => 0.5*ωn[2:end].^2);
            Cm = kron(stencil_Cm, TΔϵ);

            @info "Calculating: total matrix"
            A = blockdiag(As...) + Cp + Cm;
        else
            A = As[1];
        end

        ez = dolinearsolve(A, b, CNSym);

        @info "Processing results"

        for j = 1:nfrequencies
            ezi = ez[(j-1)*length(d.grid)+1:j*length(d.grid)];
            hxi = -1/1im/ωn[j]/μ₀*Syf[d.sharedpml ? 1 : j]*δyf*ezi;
            hyi = 1/1im/ωn[j]/μ₀*Sxf[d.sharedpml ? 1 : j]*δxf*ezi;
            fields[i, j] = FieldTM(d.grid, ωn[j], ezi, hxi, hyi)
        end
    end

    return fields
end
