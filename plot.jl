export plot_fields_abs
using PyPlot

function plot_fields_abs(mod::Modulator, Ez; cbar=false)
	epsr = real.(mod.geom.epsr)
	mod_reg = real.(mod.epsr_delta)
	xrange = mod.geom.xrange
	yrange = mod.geom.yrange
	Emx = maximum(abs.(Ez))
	N = size(epsr)
	Nfreqs = size(Ez)[1]
	fig, axes = subplots(Nfreqs,sharex=true,sharey=true,figsize=(9,12))
	for i in 1:Nfreqs
	    Z = (Ez[i,:,:].')./Emx;
	    ax=axes[i];
	    mappable = ax[:imshow](abs.(Z), cmap="inferno", extent=[xrange[1],xrange[2],yrange[1],yrange[2]]/1e-6);
	    ax[:set_ylabel](L"$y$ ($\mu$m)")
	    ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
	    	linspace(yrange[1],yrange[2],N[2])./1e-6,epsr.',1,colors=("grey",),linewidths=(.75,),alpha=1)
	    ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
	    	linspace(yrange[1],yrange[2],N[2])./1e-6,mod_reg.',1,colors=("green",),linewidths=(.75,),alpha=1)
	    if cbar
	    	colorbar(mappable,ax=ax,label=L"$\vert E \vert$")
	    end
	end
	axes[end][:set_xlabel](L"$x$ ($\mu$m)")
	show();
end

#,vmin=-maximum(abs.(Z)),vmax=maximum(abs.(Z))