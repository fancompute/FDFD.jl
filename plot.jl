import PyPlot
using PyPlot

function plot_fields_abs(Ez, eps_r, mod_reg, xrange, yrange)
	Emx = maximum(abs.(Ez))
	N = size(eps_r)
	Nfreqs = size(Ez)[1]
	fig, axes = subplots(Nfreqs,sharex=true,sharey=true,figsize=(9,12))
	for i in 1:Nfreqs
	    Z = (Ez[i,:,:].')./Emx;
	    ax=axes[i];
	    mappable = ax[:imshow](abs.(Z), cmap="inferno", extent=[xrange[1],xrange[2],yrange[1],yrange[2]]/1e-6);
	    ax[:set_ylabel](L"$y$ ($\mu$m)")
	    ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
	    	linspace(yrange[1],yrange[2],N[2])./1e-6,eps_r.',1,colors=("grey",),linewidths=(.75,),alpha=1)
	    ax[:contour](linspace(xrange[1],xrange[2],N[1])./1e-6,
	    	linspace(yrange[1],yrange[2],N[2])./1e-6,mod_reg.',1,colors=("green",),linewidths=(.75,),alpha=1)
	    colorbar(mappable,ax=ax,label=L"$\vert E \vert$")
	end
	axes[end][:set_xlabel](L"$x$ ($\mu$m)")
	show();
end

#,vmin=-maximum(abs.(Z)),vmax=maximum(abs.(Z))