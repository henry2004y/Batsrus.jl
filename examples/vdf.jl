# Velocity distribution plot in a large PIC domain.
#
#
# Hongyang Zhou, hyzhou@umich.edu 01/30/2020

using Batsrus, PyPlot, Printf, LinearAlgebra, Statistics
# For precise colorbar control
using PyCall
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
inset_axes = inset_locator.inset_axes

## Parameters
using Batsrus.UnitfulBatsrus: me, mp
const cAlfven = 253       # average Alfven speed in G8, [km/s]
const mi = 14             # average ion mass [amu]

"""
    dist_select(fnameParticle; ParticleType='e', dir=".")

Select particle in box regions centered at `xC`, `yC`, and `zC` with lengths `xL`, `yL` and
`zL`. `ParticleType` can be ['e','i'].
"""
function dist_select(
      fnameParticle, xC = -1.90, yC = 0.0, zC = -0.1, xL = 0.005, yL = 0.2, zL = 0.07;
      ParticleType = 'e', dir = ".")
   if ParticleType == 'i'
      !occursin("region0_2", fnameParticle) && @error "Check filename!"
   elseif ParticleType == 'e'
      !occursin("region0_1", fnameParticle) && @error "Check filename!"
   end

   filefield = "3d_var_region0_0_"*fnameParticle[(end - 22):end]

   nBox = 9 # number of box regions

   # Classify particles based on locations
   region = Array{Float32, 2}(undef, 6, nBox)
   region[:, 1] = [xC-xL*3/2, xC-xL/2, yC-yL/2, yC+yL/2, zC+zL/2, zC+zL*3/2]
   region[:, 2] = [xC-xL/2, xC+xL/2, yC-yL/2, yC+yL/2, zC+zL/2, zC+zL*3/2]
   region[:, 3] = [xC+xL/2, xC+xL*3/2, yC-yL/2, yC+yL/2, zC+zL/2, zC+zL*3/2]
   region[:, 4] = [xC-xL*3/2, xC-xL/2, yC-yL/2, yC+yL/2, zC-zL/2, zC+zL/2]
   region[:, 5] = [xC-xL/2, xC+xL/2, yC-yL/2, yC+yL/2, zC-zL/2, zC+zL/2]
   region[:, 6] = [xC+xL/2, xC+xL*3/2, yC-yL/2, yC+yL/2, zC-zL/2, zC+zL/2]
   region[:, 7] = [xC-xL*3/2, xC-xL/2, yC-yL/2, yC+yL/2, zC-zL*3/2, zC-zL/2]
   region[:, 8] = [xC-xL/2, xC+xL/2, yC-yL/2, yC+yL/2, zC-zL*3/2, zC-zL/2]
   region[:, 9] = [xC+xL/2, xC+xL*3/2, yC-yL/2, yC+yL/2, zC-zL*3/2, zC-zL/2]

   particle = [Array{Float32}(undef, 3, 0) for _ in 1:nBox]

   data = load(joinpath(dir, fnameParticle))

   x, y, z = eachslice(data.x, dims = 4)

   ux_ = findfirst(x->x=="ux", data.head.wname)
   uy_ = findfirst(x->x=="uy", data.head.wname)
   uz_ = findfirst(x->x=="uz", data.head.wname)

   ux = @view data.w[:, :, :, ux_]
   uy = @view data.w[:, :, :, uy_]
   uz = @view data.w[:, :, :, uz_]

   for ip in eachindex(x)
      for iR in 1:nBox
         if region[1, iR] < x[ip] < region[2, iR] &&
            region[3, iR] < y[ip] < region[4, iR] &&
            region[5, iR] < z[ip] < region[6, iR]
            particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
            break
         end
      end
   end
   region, particle
end

"""
    dist_scan(region, particle, ParticleType='i', plottype=1)

Velocity distribution plot in 9 box regions.
`plottype`: 1: uy-ux; 2: ux-uz; 3:uy-uz; 4:u⟂O-u⟂I; 5:u⟂I-u∥; 5:u⟂O-u∥
"""
function dist_scan(region, particle, ParticleType = 'i', plottype = 1; dir = ".",
      filefield::AbstractString, nbin = 60, fs = 10)
   if ParticleType == 'i'
      binRange = [[-3.0, 3.0], [-3.0, 3.0]]
   elseif ParticleType == 'e'
      binRange = [[-10.0, 12.0], [-10.0, 10.0]]
   end

   nBox = 9
   figure(figsize = (11, 6))
   for iB in 1:nBox
      if plottype ≤ 3
         ux = particle[iB][1, :] ./ cAlfven
         uy = particle[iB][2, :] ./ cAlfven
         uz = particle[iB][3, :] ./ cAlfven
      else
         dBx, dBy, dBz = GetMeanField(filefield, region[:, iB]; dir)

         dPar = [dBx; dBy; dBz] # Parallel direction
         dPerpI = cross([0; -1; 0], dPar) # Perpendicular direction in-plane
         dPerpO = cross(dPar, dPerpI) # Perpendicular direction out-of-plane

         uPar = transpose(particle[iB][1:3, :])*dPar ./ cAlfven
         uPerpI = transpose(particle[iB][1:3, :])*dPerpI ./ cAlfven
         uPerpO = transpose(particle[iB][1:3, :])*dPerpO ./ cAlfven
      end

      ax = subplot(3, 4, iB+ceil(iB/3))
      if plottype == 1
         h = hist2D(uy, ux, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype == 2
         h = hist2D(ux, uz, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype == 3
         h = hist2D(uy, uz, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype == 4
         h = hist2D(uPerpO, uPerpI, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype == 5
         h = hist2D(uPerpI, uPar, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype==6
         h = hist2D(uPerpO, uPar, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      else
         @error "Unknown PlotVType!"
      end
      grid(true)
      axis("equal")

      if iB in (1, 4, 7, 8, 9)
         if plottype == 1
            xlabel(L"u_y", fontsize = fs)
            ylabel(L"u_x", fontsize = fs)
         elseif plottype == 2
            xlabel(L"u_x", FontSize = fs)
            ylabel(L"u_z", FontSize = fs)
         elseif plottype == 3
            xlabel(L"u_y", FontSize = fs)
            ylabel(L"u_z", FontSize = fs)
         elseif plottype == 4
            xlabel(L"u_{\perp Out}", fontsize = fs)
            ylabel(L"u_{\perp In}", fontsize = fs)
         elseif plottype == 5
            xlabel(L"u_{\perp In}", FontSize = fs)
            ylabel(L"u_\parallel", FontSize = fs)
         elseif plottype == 6
            xlabel(L"u_{\perp Out}", FontSize = fs)
            ylabel(L"u_\parallel", FontSize = fs)
         end
      end
      title(@sprintf("%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]", iB, region[1, iB],
         region[2, iB], region[5, iB], region[6, iB]))
      colorbar()
      plt.set_cmap("jet")
      #clim(1e-2,10^0.3)

      if ParticleType == 'e'
         str = "electron"
      elseif ParticleType == 'i'
         str = "ion"
      end
      text(0.05, 0.05, str, FontSize = fs, transform = ax.transAxes)
   end
end

"""
    GetMeanField(filefield, limits; dir=".")

Get the average field direction in a box region specified by `limits`.
"""
function GetMeanField(filefield, limits; dir = ".")
   data = load(joinpath(dir, filefield))

   x = data.x[:, :, :, 1]
   y = data.x[:, :, :, 2]
   z = data.x[:, :, :, 3]

   bx_ = findfirst(x->x=="Bx", data.head.wname)
   by_ = findfirst(x->x=="By", data.head.wname)
   bz_ = findfirst(x->x=="Bz", data.head.wname)

   Bx = @view data.w[:, :, :, bx_]
   By = @view data.w[:, :, :, by_]
   Bz = @view data.w[:, :, :, bz_]

   xnew, ynew, znew, BxNew, ByNew, BzNew = subvolume(x, y, z, Bx, By, Bz, limits)

   # Average over the selected volume
   B̄x, B̄y, B̄z = mean(BxNew), mean(ByNew), mean(BzNew)

   # Normalize vector
   Length = √(B̄x^2 + B̄y^2 + B̄z^2)
   dBx, dBy, dBz = B̄x/Length, B̄y/Length, B̄z/Length

   dBx, dBy, dBz
end

function plotExCut(filefield::String, region, xC, yC, zC, xL, yL, zL;
      dir = ".", fs = 16, sequence = 129)
   plotrange = [xC-xL*16, xC+xL*16, zC-zL*5, zC+zL*5]
   # Sample region plot over contour
   data = load(joinpath(dir, filefield))

   bx_ = findfirst(x->x=="Bx", data.head.wname)
   bz_ = findfirst(x->x=="Bz", data.head.wname)

   Bx = @view data.w[:, :, :, bx_]
   Bz = @view data.w[:, :, :, bz_]

   subplot(3, 4, (1, 9))
   cutplot(data, "Ex"; cut = "y", sequence, plotrange)
   colorbar()
   axis("scaled")
   plt.set_cmap("RdBu_r")
   clim(-9e4, 9e4)

   streamslice(data, "Bx;Bz"; cut = "y", sequence, plotrange, color = "k", density = 1.0)

   xlabel(L"x [R_G]", fontsize = fs)
   ylabel(L"z [R_G]", fontsize = fs)
   title(L"Ex [\mu V/m]")

   nBox = 9
   for iB in 1:nBox
      rect = matplotlib.patches.Rectangle((region[1, iB], region[5, iB]),
         region[2, iB]-region[1, iB], region[6, iB]-region[5, iB],
         ec = "r", lw = 1.2, fill = false) # facecolor="none"
      ax = gca()
      ax.add_patch(rect)
   end
end

function dist_plot(pType = 'e')
   dir = "."
   nBox = 4
   nbin = 60
   fs = 10
   subplotlist = [1, 2, 3, 4, 1, 2, 3, 4]
   plottype = (1, 1, 1, 1, 3, 3, 3, 3)

   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   filefield = "3d_var_region0_0_"*fnameE[(end - 22):end]

   binRangeI = [[-3.0, 3.0], [-3.0, 3.0]]
   binRangeE = [[-7.0, 13.0], [-10.0, 10.0]]

   # Classify particles based on locations
   region = Array{Float32, 2}(undef, 6, nBox)
   particle = [Array{Float32}(undef, 3, 0) for _ in 1:nBox]

   # Electron
   if pType == 'e'
      region[:, 1] = [-1.930, -1.925, -0.08, 0.08, -0.30, -0.26]
      region[:, 2] = [-1.930, -1.925, -0.08, 0.08, -0.20, -0.16]
      region[:, 3] = [-1.930, -1.925, -0.08, 0.08, -0.10, -0.06]
      region[:, 4] = [-1.930, -1.925, -0.08, 0.08, 0.00, 0.04]

      data = load(joinpath(dir, fnameE))

      x = @view data.x[:, :, :, 1]
      y = @view data.x[:, :, :, 2]
      z = @view data.x[:, :, :, 3]

      ux_ = findfirst(x->x=="ux", data.head.wname)
      uy_ = findfirst(x->x=="uy", data.head.wname)
      uz_ = findfirst(x->x=="uz", data.head.wname)

      ux = @view data.w[:, :, :, ux_]
      uy = @view data.w[:, :, :, uy_]
      uz = @view data.w[:, :, :, uz_]

      for ip in eachindex(x)
         for iR in 1:nBox
            if region[1, iR] < x[ip] < region[2, iR] &&
               region[3, iR] < y[ip] < region[4, iR] &&
               region[5, iR] < z[ip] < region[6, iR]
               particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
               break
            end
         end
      end

      binRange = binRangeE
   elseif pType == 'i'
      region[:, 1] = [-1.930, -1.920, -0.08, 0.08, -0.30, -0.26]
      region[:, 2] = [-1.930, -1.920, -0.08, 0.08, -0.20, -0.16]
      region[:, 3] = [-1.930, -1.920, -0.08, 0.08, -0.10, -0.06]
      region[:, 4] = [-1.930, -1.920, -0.08, 0.08, 0.00, 0.04]

      data = load(joinpath(dir, fnameI))

      x = @view data.x[:, :, :, 1]
      y = @view data.x[:, :, :, 2]
      z = @view data.x[:, :, :, 3]

      ux_ = findfirst(x->x=="ux", data.head.wname)
      uy_ = findfirst(x->x=="uy", data.head.wname)
      uz_ = findfirst(x->x=="uz", data.head.wname)

      ux = @view data.w[:, :, :, ux_]
      uy = @view data.w[:, :, :, uy_]
      uz = @view data.w[:, :, :, uz_]

      for ip in eachindex(x)
         for iR in 1:nBox
            if region[1, iR] < x[ip] < region[2, iR] &&
               region[3, iR] < y[ip] < region[4, iR] &&
               region[5, iR] < z[ip] < region[6, iR]
               particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
               break
            end
         end
      end

      binRange = binRangeI
   end

   # Normalized quantities
   fig, ax = plt.subplots(4, 2, figsize = (4.2, 8.0))
   if pType == 'e'
      fig.suptitle("electron", fontsize = 14, fontweight = "bold")
   elseif pType == 'i'
      fig.suptitle("ion", fontsize = 14, fontweight = "bold")
   end

   xlPos = (0.50, -0.10)
   ylPos = (-0.25, 0.50)
   c = Vector{PyObject}(undef, nBox)
   plt.set_cmap("jet")

   for i in 1:(nBox * 2)
      ax[i].tick_params(which = "both", direction = "in", top = true, right = true)
      ax[i].minorticks_on()
      iB = subplotlist[i]
      if plottype[i] ≤ 3
         ux = particle[iB][1, :] ./ cAlfven
         uy = particle[iB][2, :] ./ cAlfven
         uz = particle[iB][3, :] ./ cAlfven
      else
         dBx, dBy, dBz = GetMeanField(filefield, region[:, iB]; dir = dir)

         dPar = [dBx; dBy; dBz]         # Parallel direction
         dPerpI = cross([0; -1; 0], dPar) # Perpendicular direction in-plane
         dPerpO = cross(dPar, dPerpI)     # Perpendicular direction out-of-plane

         uPar = transpose(particle[iB][1:3, :])*dPar ./ cAlfven
         uPerpI = transpose(particle[iB][1:3, :])*dPerpI ./ cAlfven
         uPerpO = transpose(particle[iB][1:3, :])*dPerpO ./ cAlfven
      end

      if plottype[i] == 1
         h = ax[i].hist2d(uy, ux, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype[i] == 2
         h = ax[i].hist2d(ux, uz, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype[i] == 3
         h = ax[i].hist2d(uy, uz, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype[i] == 4
         h = ax[i].hist2d(uPerpO, uPerpI, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype[i] == 5
         h = ax[i].hist2d(uPerpI, uPar, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      elseif plottype[i] == 6
         h = ax[i].hist2d(uPerpO, uPar, bins = nbin,
            norm = matplotlib.colors.LogNorm(), density = true, range = binRange)
      else
         @error "Unknown PlotVType!"
      end
      ax[i].set_aspect("equal", "box")
      #ax[i].grid(true)
      ax[i].axvline(x = 0.0, ls = "--", linewidth = 0.5, color = "k")
      ax[i].axhline(y = 0.0, ls = "--", linewidth = 0.5, color = "k")
      (i % 4 != 0) && ax[i].axes.xaxis.set_ticklabels([])

      if plottype[i] == 1
         if i in (4, 8)
            ax[i].annotate(L"u_y", xy = (0.50, -0.20), xycoords = "axes fraction",
               fontsize = fs)
         else
            ax[i].annotate(L"u_y", xy = xlPos, xycoords = "axes fraction",
               fontsize = fs)
         end
         ax[i].annotate(L"u_x", xy = ylPos, xycoords = "axes fraction",
            fontsize = fs)
      elseif plottype[i] == 2
         ax[i].annotate(L"u_x", xy = xlPos, xycoords = "axes fraction",
            fontsize = fs)
         ax[i].annotate(L"u_z", xy = ylPos, xycoords = "axes fraction",
            fontsize = fs)
      elseif plottype[i] == 3
         if i in (4, 8)
            ax[i].annotate(L"u_y", xy = (0.50, -0.20), xycoords = "axes fraction",
               fontsize = fs)
         else
            ax[i].annotate(L"u_y", xy = xlPos, xycoords = "axes fraction",
               fontsize = fs)
         end
         ax[i].annotate(L"u_z", xy = ylPos, xycoords = "axes fraction",
            fontsize = fs)
      elseif plottype[i] == 4
         ax[i].annotate(L"u_{\perp Out}", xy = xlPos, xycoords = "axes fraction",
            fontsize = fs)
         ax[i].annotate(L"u_{\perp In}", xy = ylPos, xycoords = "axes fraction",
            fontsize = fs)
      elseif plottype[i] == 5
         ax[i].annotate(L"u_{\perp In}", xy = xlPos, xycoords = "axes fraction",
            fontsize = fs)
         ax[i].annotate(L"u_\parallel", xy = ylPos, xycoords = "axes fraction",
            fontsize = fs)
      elseif plottype[i] == 6
         ax[i].annotate(L"u_{\perp Out}", xy = xlPos, xycoords = "axes fraction",
            fontsize = fs)
         ax[i].annotate(L"u_\parallel", xy = ylPos, xycoords = "axes fraction",
            fontsize = fs)
      end

      if pType == 'e'
         ax[i].text(-0.05,
            1.02,
            @sprintf("x[%3.3f,%3.3f], z[%3.3f,%3.3f]",
               region[1, iB], region[2, iB], region[5, iB], region[6, iB]),
            FontSize = 7,
            transform = ax[i].transAxes)
         ax[i].text(0.05, 0.05, string(iB), FontSize = fs, color = "r",
            transform = ax[i].transAxes)
      elseif pType == 'i'
         ax[i].text(-0.05,
            1.02,
            @sprintf("x[%3.3f,%3.3f], z[%3.3f,%3.3f]",
               region[1, iB], region[2, iB], region[5, iB], region[6, iB]),
            FontSize = 7, transform = ax[i].transAxes)
         ax[i].text(0.83, 0.05, string(iB+4), FontSize = fs, color = "c",
            transform = ax[i].transAxes)
      end
   end

   im = plt.gca().get_children()[1]
   cax = fig.add_axes([0.1, 0.93, 0.8, 0.02])
   colorbar(im, ax = ax, cax = cax, orientation = "horizontal")
   cax.tick_params(which = "both", axis = "x", direction = "out", color = "k")

   return ax
end

function show_box_region()
   dir = "/Users/hyzhou/Documents/Computer/DeepBlue"
   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   filefield = "3d_var_region0_0_"*fnameE[(end - 22):end]

   fs = 10
   plotrange = [-2.0, -1.8, -0.35, 0.35]
   sequence = 129 # cut plane index starting from -

   data = load(joinpath(dir, filefield()))

   X, Z, Bx = cutdata(data, "Bx"; dir = "y", sequence, plotrange)
   X, Z, Bz = cutdata(data, "Bz"; dir = "y", sequence, plotrange)
   X, Z, Ex = cutdata(data, "Ex"; dir = "y", sequence, plotrange)

   fig, ax = plt.subplots(1, 1, figsize = (8.0, 4.0))
   #plt.set_cmap("RdBu_r")
   plt.set_cmap("seismic")

   region = Array{Float32, 2}(undef, 6, 8)

   region[:, 1] = [-1.930, -1.925, -0.08, 0.08, -0.30, -0.26]
   region[:, 2] = [-1.930, -1.925, -0.08, 0.08, -0.20, -0.16]
   region[:, 3] = [-1.930, -1.925, -0.08, 0.08, -0.10, -0.06]
   region[:, 4] = [-1.930, -1.925, -0.08, 0.08, 0.00, 0.04]

   region[:, 5] = [-1.930, -1.920, -0.08, 0.08, -0.30, -0.26]
   region[:, 6] = [-1.930, -1.920, -0.08, 0.08, -0.20, -0.16]
   region[:, 7] = [-1.930, -1.920, -0.08, 0.08, -0.10, -0.06]
   region[:, 8] = [-1.930, -1.920, -0.08, 0.08, 0.00, 0.04]

   c = ax.contourf(Z, X, Ex ./ 1e3, 100, norm = matplotlib.colors.DivergingNorm(0),
      vmin = -12e1, vmax = 12e1)

   # create an axes on the right side of ax. The width of cax will be 5%
   # of ax and the padding between cax and ax will be fixed at 0.05 inch.
   divider = axes_grid1.make_axes_locatable(ax)
   cax = divider.append_axes("right", size = "5%", pad = 0.05)
   colorbar(c, cax = cax)
   ax.invert_yaxis()
   ax.set_aspect("equal", "box")

   ax.set_xlabel(L"z [R_G]", fontsize = 14)
   ax.set_ylabel(L"x [R_G]", fontsize = 14)
   ax.set_title(L"Ex [mV/m]", fontsize = 14)

   x = @view X[:, 1]
   z = @view Z[1, :]

   seeds = select_seeds(x[10:(end - 10)], z[10:(end - 10)]; nSeed = 1)
   xstart, zstart = seeds[1, :], seeds[2, :]
   append!(xstart, [-1.87, -1.97, -1.93, -1.93, -1.93, -1.95, -1.82, -1.93, -1.93, -1.93])
   append!(zstart, [-0.1, -0.1, -0.15, 0.0, 0.2, -0.25, -0.1, -0.08, -0.095, 0.116])

   xl = [Vector{Float32}(undef, 0) for _ in 1:length(xstart)]
   zl = [Vector{Float32}(undef, 0) for _ in 1:length(xstart)]
   for i in eachindex(xstart)
      xs, zs = xstart[i], zstart[i]
      xl[i],
      zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds = 0.02, maxstep = 20000,
         gridType = "ndgrid")
   end
   [ax.plot(zl[j], xl[j], "-", color = "k", lw = 1.0) for j in eachindex(xstart)]

   ax.contour(Z, X, Bz, [0.0], colors = "k", linestyles = "dotted", linewidths = 1.0)

   for iB in 1:4
      rect = matplotlib.patches.Rectangle((region[5, iB], region[1, iB]),
         region[6, iB]-region[5, iB], region[2, iB]-region[1, iB],
         lw = 1.0, fill = true, alpha = 0.9, fc = "r") # facecolor="none"
      ax.add_patch(rect)
   end

   for iB in 5:8
      rect = matplotlib.patches.Rectangle((region[5, iB], region[1, iB]),
         region[6, iB]-region[5, iB], region[2, iB]-region[1, iB],
         lw = 1.0, fill = true, alpha = 0.9, fc = "c") # facecolor="none"
      ax.add_patch(rect)
   end

   ax.annotate("1", (region[5, 1]-2e-2, region[2, 1]), fontsize = 10, color = "r")
   ax.annotate("2", (region[5, 2]-2e-2, region[2, 2]), fontsize = 10, color = "r")
   ax.annotate("3", (region[5, 3]-2e-2, region[2, 3]+0.5e-2), fontsize = 10, color = "r")
   ax.annotate("4", (region[5, 4]-2e-2, region[2, 4]+1e-2), fontsize = 10, color = "r")

   for i in 5:8
      ax.annotate(string(i), (region[5, i]-2e-2, region[2, i]+0.2e-2), fontsize = 10,
         color = "c")
   end

   return ax
end

function HF_velocity()
   dir = "/Users/hyzhou/"
   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   filefield = "3d_var_region0_0_"*fnameE[(end - 22):end]

   data = load(joinpath(dir, filefield))

   x = data.x[:, :, :, 1]
   y = data.x[:, :, :, 2]
   z = data.x[:, :, :, 3]

   ux_ = findfirst(x->lowercase(x)=="uxs1", data.head.wname)
   uy_ = findfirst(x->lowercase(x)=="uys1", data.head.wname)
   uz_ = findfirst(x->lowercase(x)=="uzs1", data.head.wname)

   Ux = @view data.w[:, :, :, ux_]
   Uy = @view data.w[:, :, :, uy_]
   Uz = @view data.w[:, :, :, uz_]

   limits = [-1.910, -1.900, -0.08, 0.08, 0.11, 0.15]
   xnew, ynew, znew, UxNew, UyNew, UzNew = subvolume(x, y, z, Ux, Uy, Uz, limits)

   # Average over the selected volume
   Ūx, Ūy, Ūz = mean(UxNew), mean(UyNew), mean(UzNew)

   return Ūx, Ūy, Ūz
end

dir = "/Users/hyzhou/Documents/Computer/DeepBlue"
filefield = "3d_var_region0_0_t00001640_n00020369.out"
PType = 'e'
plottype = 2

if PType == 'e'
   fnameParticle = "cut_particles0_region0_1_t00001640_n00020369.out"
elseif PType == 'i'
   fnameParticle = "cut_particles1_region0_2_t00001640_n00020369.out"
end

# Define regions
xC, yC, zC = -1.75, 0.0, -0.2
xL, yL, zL = 0.008, 0.2, 0.03 # box length in x,y,z

@time region,
particle = dist_select(
   fnameParticle, xC, yC, zC, xL, yL, zL,
   dir = dir, ParticleType = PType)

@time dist_scan(region, particle, PType, plottype; dir, filefield)

@time plotExCut(filefield, region, xC, yC, zC, xL, yL, zL; dir)

#ax = dist_plot('i')
#ux, uy, uz = HF_velocity()
#show_box_region() # horizontal
