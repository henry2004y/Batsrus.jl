# Script for plotting a series of cuts.
#
# Hongyang Zhou, hyzhou@umich.edu

using Batsrus, PyPlot, Glob

matplotlib.rc("image", cmap=:turbo)

filenames = glob("GM/IO2/z*out")

dir = "GM/IO2"
var = "rho"
levels = 50
plotinterval = 0.05
plotrange = [-Inf,Inf,-Inf,Inf]
vmin = 0.0
vmax = 10.0
colorscale = :linear

fig, ax = subplots(figsize=(10,7.5))

data = readdata(basename(filenames[1]); dir)
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

nfiles = length(filenames)
if nfiles > 1
   for i in 2:nfiles
      f = filenames[i] |> basename
      @info "$i / $nfiles, $f"
      local data = readdata(f; dir)

      contourf(data, var, levels; ax, plotrange, plotinterval, innermask=true, norm=cnorm)
   
      title("2D Magnetosphere, MHD, t=$(data.head.time)s")
      xlabel("X [Re]", fontsize=14)
      ylabel("Y [Re]", fontsize=14)

      plt.savefig("out/"*lpad(i, 4, '0')*".png", bbox_inches="tight")
      plt.cla()
   end
end

plt.close()