"    probe_field(field::Field, xy::AbstractArray)"
function probe_field(field::Field, xy::AbstractArray)
    (indx, indy) = coord2ind(field.grid, xy)
    isa(field, FieldTM) && return field.Ez[indx, indy]
    isa(field, FieldTE) && return field.Hz[indx, indy]
end

"    scattering_parameters(fields::Array{<:Field}, d::AbstractDevice)"
function scattering_parameters(fields::Array{<:Field}, d::AbstractDevice)
    @assert length(fields) == length(d.ω)
    (ϵ₀, μ₀, c₀) = normalize_parameters(d)
    S = zeros(Complex, length(d.ω), length(d.modes))
    for i = 1:length(d.ω)
        ω = d.ω[i];
        for j = 1:length(d.modes)
            mode = d.modes[j];
            (β, Ej, indx, indy) = get_modes(d, mode.pol, ω, mode.neff, 1, mode.coor, mode.dir, mode.width)
            Ec = fields[i].Ez[indx,indy] # computed field
            Z = β/ω/ϵ₀

            if j == 1
                numer = (Ec-Ej).*conj.(Ej)
            else
                numer = Ec.*conj.(Ej)
            end
            denom = Ej.*conj.(Ej)
            S[i, j] = sum(numer)/sum(denom)
        end
    end
    return S
end

"    poynting(field::Field)"
function poynting(field::Field)
    if isa(field, FieldTM)
        Ez_x = grid_average(field.Ez, DirectionX);
        Ez_y = grid_average(field.Ez, DirectionY);
        Sx = -0.5*real.(Ez_x.*conj(field.Hy));
        Sy =  0.5*real.(Ez_y.*conj(field.Hx));
        return Flux2D(field.grid, Sx, Sy)
    end
    if isa(field, FieldTE)
        Hz_x = grid_average(field.Hz, DirectionX);
        Hz_y = grid_average(field.Hz, DirectionY);
        Sx =  0.5*real.(field.Ey.*conj(field.Hz_x));
        Sy = -0.5*real.(field.Ex.*conj(field.Hz_y));
        return Flux2D(field.grid, Sx, Sy)
    end
    print_warn("poynting() detected neither a TM nor a TE field polarization");
end

"    flux_surface(fields::Array{<:Field}, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction)"
function flux_surface(fields::Array{<:Field}, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction)
    Nfields = length(fields);
    P = zeros(Float, Nfields);
    for i in eachindex(fields)
        P[i] = flux_surface(poynting(fields[i]), ptmid, width, nrm);
    end
    return P
end

"   flux_surface(field::Field, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction) = flux_surface(poynting(field), ptmid, width, nrm)"
flux_surface(field::Field, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction) = flux_surface(poynting(field), ptmid, width, nrm);

"    flux_surface(flux::Flux2D, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction)"
function flux_surface(flux::Flux2D, ptmid::AbstractArray{<:Real}, width::Real, nrm::Direction)
    (x0, y0) = coord2ind(flux.grid, ptmid);
    if nrm == DirectionX
        indx = x0;
        if isinf(width)
            indy = 1:flux.grid.N[2];
        else
            y1 = y2ind(flux.grid,ptmid[2]-width/2);
            y2 = y2ind(flux.grid,ptmid[2]+width/2);
            indy = y1:1:y2;
        end
        return sum(flux.Sx[indx,indy])*dy(flux.grid)
    elseif nrm == DirectionY
        error("Not implemented yet...")
    end
end