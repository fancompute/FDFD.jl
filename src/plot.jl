export plot_field

using PyPlot

function plot_field(field::Field; cbar::Bool=false, funcz=real)
	isa(field, FieldTM) && (Z = funcz.(field.Ez).');
	isa(field, FieldTE) && (Z = funcz.(field.Hz).');

	if funcz == abs
		vmin = 0;
		vmax = +maximum(abs.(Z));
		cmap = "inferno";
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
	
	fig, ax = subplots(1, sharex=true, sharey=true)
	mappable = ax[:imshow](Z, cmap=cmap, extent=extents, origin="lower", vmin=vmin, vmax=vmax);
	cbar && colorbar(mappable, ax=ax, label=L"$\vert E \vert$");
	ax[:set_xlabel](L"$x$");
	ax[:set_ylabel](L"$y$");
	show();
end