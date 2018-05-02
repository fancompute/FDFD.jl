export plot_field, plot_device, add_scalebar

using PyPlot, PyCall

"    plot_field(fields::Array{<:Field}; cbar::Bool=false, funcz=real)"
function plot_field(fields::Array{<:Field}; cbar::Bool=false, funcz=real)
	nfields = length(fields);
	fig, axs = subplots(nfields, 1, constrained_layout=true);
	for i in eachindex(fields)
		plot_field(axs[i], fields[i]; cbar=cbar, funcz=funcz);
	end
	return axs
end

"    plot_field(field::Field; cbar::Bool=false, funcz=real)"
function plot_field(field::Field; cbar::Bool=false, funcz=real)
	fig, ax = subplots(1);
	plot_field(ax, field; cbar=cbar, funcz=funcz);
	return ax
end

"    plot_field(ax::PyObject, field::Field; cbar::Bool=false, funcz=real)"
function plot_field(ax::PyObject, field::Field; cbar::Bool=false, funcz=real)
	isa(field, FieldTM) && (Z = funcz.(field.Ez).');
	isa(field, FieldTE) && (Z = funcz.(field.Hz).');

	if funcz == abs
		vmin = 0;
		vmax = +maximum(abs.(Z));
		cmap = "magma";
	elseif funcz == real
		Zmx  = maximum(abs.(Z));
		vmin = -maximum(Zmx);
		vmax = +maximum(Zmx);
		cmap = "RdBu";
	else
		error("Unknown function specified.");
	end

	extents = [ field.grid.bounds[1][1], field.grid.bounds[2][1], 
	            field.grid.bounds[1][2], field.grid.bounds[2][2] ];
	
	mappable = ax[:imshow](Z, cmap=cmap, extent=extents, origin="lower", vmin=vmin, vmax=vmax);
	cbar && colorbar(mappable, ax=ax, label=L"$\vert E \vert$");
	ax[:set_xlabel](L"$x$");
	ax[:set_ylabel](L"$y$");
end

"    plot_device(device::AbstractDevice; outline::Bool=false)"
function plot_device(device::AbstractDevice; outline::Bool=false)
	fig, ax = subplots(1);
	plot_device(ax, device; outline=outline);
	return ax
end

"    plot_device(ax::Array{<:PyObject}, device::AbstractDevice; outline::Bool=false, lc::String=\"k\", lcm::String=\"k\")"
function plot_device(axs::Array{<:PyObject}, device::AbstractDevice; outline::Bool=false, lc::String="k", lcm::String="k")
	for i in eachindex(axs)
		plot_device(axs[i], device, outline=outline, lc=lc, lcm=lcm);
	end
end

"    plot_device(ax::PyObject, device::AbstractDevice; outline::Bool=false, lc::String=\"k\", lcm::String=\"k\")"
function plot_device(ax::PyObject, device::AbstractDevice; outline::Bool=false, lc::String="k", lcm::String="k")
	Z = real.(device.ϵᵣ)';

	if outline
		ax[:contour](xc(device.grid), yc(device.grid), Z, levels=1, linewidths=0.5, colors=lc);
		axis("image");
	else
		extents = [ device.grid.bounds[1][1], device.grid.bounds[2][1], 
		            device.grid.bounds[1][2], device.grid.bounds[2][2] ];
		ax[:imshow](Z, extent=extents, origin="lower", cmap="YlGnBu");
	end

	if isa(device, ModulatedDevice)
		Z2 = abs.(device.Δϵᵣ)';
		ax[:contour](xc(device.grid), yc(device.grid), Z2, levels=1, linewidths=0.5, colors=lcm);
	end
	ax[:set_xlabel](L"$x$");
	ax[:set_ylabel](L"$y$");
end

"    add_scalebar(axs::Array{<:PyObject}, xy::AbstractArray; fc::String=\"k\", width::Number=1.0, height::Number=0.125, hide::Bool=true)"
function add_scalebar(axs::Array{<:PyObject}, xy::AbstractArray; fc::String="k", width::Number=1.0, height::Number=0.125, hide::Bool=true)
	for i in eachindex(axs)
		add_scalebar(axs[i], xy, fc=fc, width=width, height=height, hide=hide);
	end
end

"    add_scalebar(ax::PyObject, xy::AbstractArray; fc::String=\"k\", width::Number=1.0, height::Number=0.125, hide::Bool=true)"
function add_scalebar(ax::PyObject, xy::AbstractArray; fc::String="k", width::Number=1.0, height::Number=0.125, hide::Bool=true)
    ax[:add_patch](matplotlib[:patches][:Rectangle](xy-[width/2,height/2],width,height,fc=fc));
    ax[:annotate](xy=xy+[0,height/2],s=@sprintf("%d",width)*L" $\mu$m",fontsize="smaller",ha="center",va="bottom",multialignment="center",color=fc);
    if hide
    	ax[:axes][:get_xaxis]()[:set_visible](false);
    	ax[:axes][:get_yaxis]()[:set_visible](false);
    end
end