using AxisArrays, IntervalSets

export probe_field, poynting, flux_surface_integral

"     probe_field(field::Field, component::Symbol, pt:Point)"
function probe_field(field::Field, component::Symbol, pt::Point)
    xind = AxisArrays.axisindexes(axes(field)[XX], atvalue(pt.x, atol=dx(field.grid)/2))
    yind = AxisArrays.axisindexes(axes(field)[YY], atvalue(pt.y, atol=dy(field.grid)/2))
    isa(field, FieldTM) && ~in(component, [:Ez, :Hx, :Hy]) && error("$component is invalid for the TM polarization. Valid options are :Ez, :Hx, :Hy")
    isa(field, FieldTE) && ~in(component, [:Hz, :Ex, :Ey]) && error("$component is invalid for the TE polarization. Valid options are :Hz, :Ex, :Ey")
    return field[xind, yind, component]
end

"    poynting(field::Field)"
function poynting(field::Field)
    if isa(field, FieldTM)
        Ez_x = grid_average(field[:, :, :Ez], x̂);
        Ez_y = grid_average(field[:, :, :Ez], ŷ);
        Sx = -0.5*real.(Ez_x.*conj(field[:, :, :Hy]));
        Sy =  0.5*real.(Ez_y.*conj(field[:, :, :Hx]));
        return Flux(field.grid, field.ω, Sx, Sy)
    end
    if isa(field, FieldTE)
        Hz_x = grid_average(field[:, :, :Hz], x̂);
        Hz_y = grid_average(field[:, :, :Hz], ŷ);
        Sx =  0.5*real.(field[:, :, :Ey].*conj(field.Hz_x));
        Sy = -0.5*real.(field[:, :, :Ex].*conj(field.Hz_y));
        return Flux(field.grid, field.ω, Sx, Sy)
    end
end

""""
    flux_surface_integral(field::Field, center::Point, width::Real, nrm::Direction)
Compute integral of flux normal to a surface defined by `center`, `width`, and
 `normal`.
"""
function flux_surface_integral(field::Field, center::Point, width::Real, normal::Direction)
    if normal == x̂
        if isa(field, FieldTM)
            xind = AxisArrays.axisindexes(AxisArrays.axes(field)[XX], atvalue(center.x, atol=dx(field.grid)/2))
            yint = (center.y-width)..(center.y+width)
            yind = AxisArrays.axisindexes(AxisArrays.axes(field)[YY], yint)
            Ez   = view(field.data, xind:xind+1, yind, :Ez)
            Ez_x̂ = grid_average(Ez, x̂)[1,:]
            Hy   = view(field.data, xind, yind, :Hy);
            return sum(-0.5*real.(Ez_x̂.*conj(Hy)))*dy(field.grid)
        end
        if isa(field, FieldTE)
            xind = AxisArrays.axisindexes(AxisArrays.axes(field)[XX], atvalue(pt.x, atol=dx(field.grid)/2))
            yint = (center.y-width)..(center.y+width)
            yind = AxisArrays.axisindexes(AxisArrays.axes(field)[YY], yint)
            Hz   = view(field.data, xind:xind+1, yind, :Hz)
            Hz_x̂ = grid_average(Hz, x̂)[1,:]
            Ey   = view(field.data, xind, yind, :Hy);
            return sum(-0.5*real.(Ey.*conj(Hz_x̂)))*dy(field.grid)
        end
    end
    if normal == ŷ
        if isa(field, FieldTM)
            error("ŷ normal TM calculation not yet implemented")
        end
        if isa(field, FieldTE)
            error("ŷ normal TE calculation not yet implemented")
        end
    end
end

# "    scattering_parameters(fields::Array{<:Field}, d::AbstractDevice)"
# function scattering_parameters(fields::Array{<:Field}, d::AbstractDevice)
#     @assert length(fields) == length(d.ω)
#     (ϵ₀, μ₀, c₀) = normalize_parameters(d)
#     S = zeros(Complex, length(d.ω), length(d.modes))
#     for i = 1:length(d.ω)
#         ω = d.ω[i];
#         for j = 1:length(d.modes)
#             mode = d.modes[j];
#             (β, Ej, indx, indy) = get_modes(d, mode.pol, ω, mode.neff, 1, mode.coor, mode.dir, mode.width)
#             Ec = fields[i].Ez[indx,indy] # computed field
#             Z = β/ω/ϵ₀
#
#             if j == 1
#                 numer = (Ec-Ej).*conj.(Ej)
#             else
#                 numer = Ec.*conj.(Ej)
#             end
#             denom = Ej.*conj.(Ej)
#             S[i, j] = sum(numer)/sum(denom)
#         end
#     end
#     return S
# end
