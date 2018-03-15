export Geometry, Modulator
export M, dx, dy, xc, yc, xe, ye, coord2ind

mutable struct Geometry1D
	N::Integer
	xrange::Tuple{Real,Real}
	epsr::Array{Complex,1}

	function Geometry1D(xrange::Tuple{Real,Real}, epsr::Array{Complex,1})
		return new(length(epsr), xrange, epsr)
	end
end

mutable struct Geometry
	N::Tuple{Integer,Integer}
	Npml::Tuple{Integer,Integer}
	xrange::Tuple{Real,Real}
	yrange::Tuple{Real,Real}
	epsr::Array{Complex,2}
	src::Array{Complex,2}

	function Geometry(N::Tuple{Integer,Integer}, Npml::Tuple{Integer,Integer}, xrange::Tuple{Real,Real}, yrange::Tuple{Real,Real})
		epsr = ones(Complex128, N)
		src = zeros(Complex128, N)
		return new(N, Npml, xrange, yrange, epsr, src)
	end

	function Geometry(N::Tuple{Integer,Integer}, xrange::Tuple{Real,Real}, yrange::Tuple{Real,Real})
		return Geometry(N, (0, 0), xrange, yrange)
	end

	function Geometry(dh::Real, Npml::Tuple{Integer,Integer}, xrange::Tuple{Real,Real}, yrange::Tuple{Real,Real})
		Nx = Int64(round((xrange[2]-xrange[1])/dh));
		Ny = Int64(round((yrange[2]-yrange[1])/dh));
		N = (Nx, Ny);
		println(" # Generated grid size: ", N);
		return Geometry(N, Npml, xrange, yrange);
	end

	function Geometry(dh::Real, xrange::Tuple{Real,Real}, yrange::Tuple{Real,Real})
		return Geometry(dh, (0, 0), xrange, yrange);
	end	
end

mutable struct Modulator
	geom::Geometry
	Omega::Real
	Nsb::Integer
	epsr_delta::Array{Complex,2}
	epsr_delta_phi::Array{Real,2}

	function Modulator(geom::Geometry, Omega::Real, Nsb::Integer)
		epsr_delta = zeros(Complex128, geom.N)
		epsr_delta_phi = zeros(Float64, geom.N)
		return new(geom, Omega, Nsb, epsr_delta, epsr_delta_phi)
	end
end

function M(geom::Geometry)
	return prod(geom.N)
end

function dx(geom::Geometry1D)
	return (geom.xrange[2]-geom.xrange[1])/geom.N[1]
end

function dx(geom::Geometry)
	return (geom.xrange[2]-geom.xrange[1])/geom.N[1]
end

function dy(geom::Geometry)
	return (geom.yrange[2]-geom.yrange[1])/geom.N[2]
end

function xc(geom::Geometry)
	return geom.xrange[1]+dx(geom)*(0.5:1:geom.N[1])
end

function yc(geom::Geometry)
	return geom.yrange[1]+dy(geom)*(0.5:1:geom.N[2])
end

function xe(geom::Geometry)
	return geom.xrange[1]+dx(geom)*(0:1:geom.N[1])
end

function ye(geom::Geometry)
	return geom.yrange[1]+dy(geom)*(0:1:geom.N[2])
end

function coord2ind(geom::Geometry, xy)
    indx = Int64(round((xy[1]-geom.xrange[1])/(geom.xrange[2]-geom.xrange[1])*geom.N[1])+1); 
    indy = Int64(round((xy[2]-geom.yrange[1])/(geom.yrange[2]-geom.yrange[1])*geom.N[2])+1);
    return (indx, indy)
end