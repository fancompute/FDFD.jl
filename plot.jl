export plot_fields, plot_fields_abs, plot_src, plot_ϵᵣ, plot_Δϵᵣ, plot_Δϵᵣ_angle
using PyPlot

function plot_fields(geom::Geometry, Ez; cbar=false, reabs = "abs")
	ϵᵣ = real.(geom.ϵᵣ);
	xrange = geom.xrange;
	yrange = geom.yrange;
	N = size(ϵᵣ);

	fig, ax = subplots(1, sharex=true, sharey=true)

	if reabs == "abs"
		Z = abs.(Ez).';
		vmin = 0;
		vmax = +maximum(abs.(Z));
		cmap = "inferno";
		lc = "grey";
	elseif reabs == "re"
		Z = real.(Ez).';
		Zmx = maximum(abs.(Ez));
		vmin = -maximum(Zmx);
		vmax = +maximum(Zmx);
		cmap = "RdBu";
		lc = "black";
	end

	mappable = ax[:imshow](Z, cmap=cmap, extent=[xrange[1],xrange[2],yrange[1],yrange[2]]/1e-6, origin="lower", vmin=vmin, vmax=vmax);
	
	ax[:set_ylabel](L"$y$ ($\mu$m)");
	ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
		linspace(yrange[1],yrange[2],N[2])./1e-6,ϵᵣ.',1,colors=(lc,),linewidths=(.75,),alpha=1);
	if cbar
	 	colorbar(mappable, ax=ax, label=L"$\vert E \vert$");
	end
	ax[:set_xlabel](L"$x$ ($\mu$m)");
	show();
end


function plot_fields_abs(mod::Modulator, Ez; cbar=false)
	ϵᵣ = real.(mod.geom.ϵᵣ);
	Δϵᵣ = abs.(mod.Δϵᵣ);
	xrange = mod.geom.xrange;
	yrange = mod.geom.yrange;
	Emx = maximum(abs.(Ez));
	N = size(ϵᵣ);
	Nfreqs = size(Ez)[1]
	fig, axes = subplots(Nfreqs, sharex=true, sharey=true)
	for i in 1:Nfreqs
	    Z = (Ez[i,:,:].')./Emx;
	    ax=axes[i];
	    mappable = ax[:imshow](abs.(Z), cmap="inferno", extent=[xrange[1],xrange[2],yrange[1],yrange[2]]/1e-6, origin="lower");
	    ax[:set_ylabel](L"$y$ ($\mu$m)");
	    ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
	    	linspace(yrange[1],yrange[2],N[2])./1e-6,ϵᵣ.',1,colors=("grey",),linewidths=(.75,),alpha=1);
	    ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
	    	linspace(yrange[1],yrange[2],N[2])./1e-6,Δϵᵣ.',1,colors=("green",),linewidths=(.75,),alpha=1);
	    if cbar
	    	colorbar(mappable, ax=ax, label=L"$\vert E \vert$");
	    end
	end
	axes[end][:set_xlabel](L"$x$ ($\mu$m)");
	show();
end

function plot_src(geom::Geometry2D)
	xy = [geom.xrange[1], geom.xrange[2], geom.yrange[1], geom.yrange[2]]
	plt_do_color(xc(geom), yc(geom), abs.(geom.src));
	plt_do_outline(xc(geom), yc(geom), real.(geom.ϵᵣ));
end

function plot_ϵᵣ(geom::Geometry2D)
	xy = [geom.xrange[1], geom.xrange[2], geom.yrange[1], geom.yrange[2]];
	plt_do_color(xc(geom), yc(geom), real.(geom.ϵᵣ));
end

function plot_Δϵᵣ(mod::Modulator)
	geom = mod.geom;
	plt_do_color(xc(geom), yc(geom), abs.(mod.Δϵᵣ));
end

function plot_Δϵᵣ_angle(mod::Modulator)
	geom = mod.geom;
	plt_do_color(xc(geom), yc(geom), angle.(mod.Δϵᵣ)/pi);
end

function plt_do_color(x, y, quantity)
	pcolormesh(x/1e-6, y/1e-6, quantity.');
	axis("scaled");
	ylabel(L"$y$ ($\mu$m)");
	xlabel(L"$x$ ($\mu$m)");
end

function plt_do_outline(x, y, quantity; lc="grey")
	contour(x/1e-6, y/1e-6, quantity.',levels=[maximum(quantity)],colors=(lc,),linewidths=(.75,),alpha=1);
	axis("scaled");
	ylabel(L"$y$ ($\mu$m)");
	xlabel(L"$x$ ($\mu$m)");
end