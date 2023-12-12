# Script for plotting a series of cuts.
#
# Hongyang Zhou, hyzhou@umich.edu

using Batsrus, PyPlot, Glob

dir = "RESULTS/run_2ndOrder_10s_PIC_10x10_xy2d/GM"
var = "rho"
plotinterval = 0.05
plotrange = [-Inf,Inf,-Inf,Inf]
vmin = 0.0
vmax = 10.0
colorscale = :linear
levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)

filenames = glob(joinpath(dir,"z*out"))
nfiles = length(filenames)

fig, ax = subplots(figsize=(10,7.5))

@info "1 / $nfiles, $(basename(filenames[1]))"
data = load(filenames[1])
varIndex_ = Batsrus.findindex(data, var)
cnorm, cticks = @views Batsrus.set_colorbar(colorscale, vmin, vmax, data.w[:,:,varIndex_])

c = contourf(data, var, levels; ax, plotrange, plotinterval, innermask=true, norm=cnorm)
title("2D Magnetosphere, MHD, t=$(data.head.time)s")

axis("scaled")
xlabel("X [Re]", fontsize=14)
ylabel("Y [Re]", fontsize=14)

cb = colorbar(c; ax, ticks=cticks)
cb_title = cb.ax.set_ylabel("Density, [amu/cc]", fontsize=14)
plt.savefig("out/"*lpad(1, 4, '0')*".png", bbox_inches="tight")
plt.cla()

if nfiles > 1
   for i in 2:nfiles
      f = filenames[i]
      @info "$i / $nfiles, $f"
      local data = load(f)

      contourf(data, var, levels; ax, plotrange, plotinterval, innermask=true, norm=cnorm)
   
      title("2D Magnetosphere, MHD, t=$(round(data.head.time, digits=1))s")
      xlabel("X [Re]", fontsize=14)
      ylabel("Y [Re]", fontsize=14)

      plt.savefig("out/"*lpad(i, 4, '0')*".png", bbox_inches="tight")
      plt.cla()
   end
end

plt.close()