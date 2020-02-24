# Distribution plot in a large PIC domain.
#
#
# Hongyang Zhou, hyzhou@umich.edu 01/30/2020

using VisAna, PyPlot, Printf, LinearAlgebra, Statistics
# For precise colorbar control
using PyCall
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
inset_axes = inset_locator.inset_axes

## Parameters
const cAlfven = 253       # average Alfven speed in G8, [km/s]
const me = 9.10938356e-31 # electron mass, [kg]
const mp = 1.6726219e-27  # proton mass, [kg]
const mi = 14             # average ion mass [amu]
#const nBox = 9            # number of box regions

"""
	dist_select(fnameParticle; ParticleType='e', dir=".")

Select particle in regions.
`ParticleType` in ['e','i'].
xC, yC, zC are the central box center position.
"""
function dist_select(fnameParticle, xC=-1.90, yC=0.0, zC=-0.1,
   xL=0.005, yL=0.2, zL=0.07; ParticleType='e',
   dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia/",)

   if ParticleType == 'i'
      !occursin("region0_2", fnameParticle) && @error "Check filename!"
   elseif ParticleType == 'e'
      !occursin("region0_1", fnameParticle) && @error "Check filename!"
   end

   fnameField = "3d_var_region0_0_"*fnameParticle[end-22:end]

   nBox = 9

   # Classify particles based on locations
   region = Array{Float32,2}(undef,6,nBox)
   region[:,1] = [xC-xL*3/2, xC-xL/2,   yC-yL/2, yC+yL/2, zC+zL/2,   zC+zL*3/2]
   region[:,2] = [xC-xL/2,   xC+xL/2,   yC-yL/2, yC+yL/2, zC+zL/2,   zC+zL*3/2]
   region[:,3] = [xC+xL/2,   xC+xL*3/2, yC-yL/2, yC+yL/2, zC+zL/2,   zC+zL*3/2]
   region[:,4] = [xC-xL*3/2, xC-xL/2,   yC-yL/2, yC+yL/2, zC-zL/2,   zC+zL/2]
   region[:,5] = [xC-xL/2,   xC+xL/2,   yC-yL/2, yC+yL/2, zC-zL/2,   zC+zL/2]
   region[:,6] = [xC+xL/2,   xC+xL*3/2, yC-yL/2, yC+yL/2, zC-zL/2,   zC+zL/2]
   region[:,7] = [xC-xL*3/2, xC-xL/2,   yC-yL/2, yC+yL/2, zC-zL*3/2, zC-zL/2]
   region[:,8] = [xC-xL/2,   xC+xL/2,   yC-yL/2, yC+yL/2, zC-zL*3/2, zC-zL/2]
   region[:,9] = [xC+xL/2,   xC+xL*3/2, yC-yL/2, yC+yL/2, zC-zL*3/2, zC-zL/2]

   particle = [Array{Float32}(undef, 3, 0) for _ in 1:nBox]

   head, data = readdata(fnameParticle, dir=dir)

   x = @view data[1].x[:,:,:,1]
   y = @view data[1].x[:,:,:,2]
   z = @view data[1].x[:,:,:,3]

   ux_ = findfirst(x->x=="ux", head[1][:wnames])
   uy_ = findfirst(x->x=="uy", head[1][:wnames])
   uz_ = findfirst(x->x=="uz", head[1][:wnames])

   ux = @view data[1].w[:,:,:,ux_]
   uy = @view data[1].w[:,:,:,uy_]
   uz = @view data[1].w[:,:,:,uz_]

   for ip = 1:length(x)
      for iR = 1:nBox
         if region[1,iR] < x[ip] < region[2,iR] &&
            region[3,iR] < y[ip] < region[4,iR] &&
            region[5,iR] < z[ip] < region[6,iR]

            particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
            break
         end
      end
   end
   return region, particle
end

"""
	dist_scan(region, particle, ParticleType='i', PlotVType=1)

Velocity distribution plot in 9 regions.
`PlotVType`: 1: uy-ux; 2: ux-uz; 3:uy-uz; 4:u⟂O-u⟂I; 5:u⟂I-u∥; 5:u⟂O-u∥
"""
function dist_scan(region, particle, ParticleType='i', PlotVType=1; dir=".",
   fnameField::String, nbin=60, fs=10)

   if ParticleType == 'i'
      binRange = [[-3.,3.], [-3.,3.]]
   elseif ParticleType == 'e'
      binRange = [[-10.,12.], [-10.,10.]]
   end

   nBox = 9
   figure(figsize=(11,6))
   for iB = 1:nBox
      if PlotVType ≤ 3
         ux = particle[iB][1,:] ./ cAlfven
         uy = particle[iB][2,:] ./ cAlfven
         uz = particle[iB][3,:] ./ cAlfven
      else
         dBx, dBy, dBz = GetMeanField(fnameField, region[:,iB]; dir=dir)

         dPar = [dBx; dBy; dBz] # Parallel direction
         dPerpI = cross([0; -1; 0], dPar) # Perpendicular direction in-plane
         dPerpO = cross(dPar, dPerpI) # Perpendicular direction out-of-plane

         uPar = transpose(particle[iB][1:3,:])*dPar ./ cAlfven
         uPerpI = transpose(particle[iB][1:3,:])*dPerpI ./ cAlfven
         uPerpO = transpose(particle[iB][1:3,:])*dPerpO ./ cAlfven
      end

      ax = subplot(3,4,iB+ceil(iB/3))
      if PlotVType==1
         h = hist2D(uy, ux, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==2
         h = hist2D(ux, uz, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==3
         h = hist2D(uy, uz, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==4
         h = hist2D(uPerpO, uPerpI, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==5
         h = hist2D(uPerpI, uPar, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==6
         h = hist2D(uPerpO, uPar, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      else
         @error "Unknown PlotVType!"
      end
      grid(true)
      axis("equal")

      if iB in (1,4,7,8,9)
         if PlotVType==1
            xlabel(L"u_y",fontsize=fs)
            ylabel(L"u_x",fontsize=fs)
         elseif PlotVType==2
            xlabel(L"u_x",FontSize=fs)
            ylabel(L"u_z",FontSize=fs)
         elseif PlotVType==3
            xlabel(L"u_y",FontSize=fs)
            ylabel(L"u_z",FontSize=fs)
         elseif PlotVType==4
            xlabel(L"u_{\perp Out}",fontsize=fs)
            ylabel(L"u_{\perp In}",fontsize=fs)
         elseif PlotVType==5
            xlabel(L"u_{\perp In}",FontSize=fs)
            ylabel(L"u_\parallel",FontSize=fs)
         elseif PlotVType==6
            xlabel(L"u_{\perp Out}",FontSize=fs)
            ylabel(L"u_\parallel",FontSize=fs)
         end
      end
      title(@sprintf("%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]",iB,region[1,iB],
         region[2,iB],region[5,iB],region[6,iB]))
      colorbar()
      plt.set_cmap("jet")
      #clim(1e-2,10^0.3)

      if ParticleType == 'e'
         str = "electron"
      elseif ParticleType == 'i'
         str = "ion"
      end
      text(0.05,0.05,str, FontSize=fs, transform=ax.transAxes)
   end

end

"""
	GetMeanField(fnameField, limits; dir=".")

GetMeanField Get the average field direction in limited region.
   * Extract the average field from field data
"""
function GetMeanField(fnameField::String, limits; dir=".")

   # Get the average field direction in limited region
   head, data = readdata(fnameField, dir=dir)

   x = data[1].x[:,:,:,1]
   y = data[1].x[:,:,:,2]
   z = data[1].x[:,:,:,3]

   bx_ = findfirst(x->x=="Bx", head[1][:wnames])
   by_ = findfirst(x->x=="By", head[1][:wnames])
   bz_ = findfirst(x->x=="Bz", head[1][:wnames])

   Bx = @view data[1].w[:,:,:,bx_]
   By = @view data[1].w[:,:,:,by_]
   Bz = @view data[1].w[:,:,:,bz_]

   xnew, ynew, znew, BxNew, ByNew, BzNew = subvolume(x,y,z, Bx,By,Bz, limits)

   # Average over the selected volume
   B̄x, B̄y, B̄z = mean(BxNew), mean(ByNew), mean(BzNew)

   # Normalize vector
   Length = √(B̄x^2 + B̄y^2 + B̄z^2)
   dBx, dBy, dBz = B̄x/Length, B̄y/Length, B̄z/Length

   return dBx, dBy, dBz
end



function plotExCut(fnameField::String, region, xC, yC, zC, xL, yL, zL;
   dir="/Users/hyzhou", fs=16, cutPlane=129)

   plotrange = [xC-xL*16, xC+xL*16, zC-zL*5, zC+zL*5]
   # Sample region plot over contour
   head, data = readdata(fnameField, dir=dir)

   bx_ = findfirst(x->x=="Bx", head[1][:wnames])
   bz_ = findfirst(x->x=="Bz", head[1][:wnames])

   Bx = @view data[1].w[:,:,:,bx_]
   Bz = @view data[1].w[:,:,:,bz_]

   subplot(3,4,(1,9))
   cutplot(data[1],head[1],"Ex",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)
   colorbar()
   axis("scaled")
   plt.set_cmap("RdBu_r")
   clim(-9e4,9e4)

   streamslice(data[1],head[1],"Bx;Bz",cut='y',cutPlaneIndex=cutPlane,
      color="k", density=1.0, plotrange=plotrange)

   xlabel(L"x [R_G]", fontsize=fs)
   ylabel(L"z [R_G]", fontsize=fs)
   title(L"Ex [\mu V/m]")

   nBox = 9
   for iB = 1:nBox
      rect = matplotlib.patches.Rectangle( (region[1,iB], region[5,iB]),
      region[2,iB]-region[1,iB], region[6,iB]-region[5,iB],
      ec="r", lw=1.2, fill=false) # facecolor="none"
      ax = gca()
      ax.add_patch(rect)
   end

   #=
   # streamline function requires the meshgrid format strictly
   s = streamslice(cut1",cut2",Bx",Bz",1,"linear")
   for is = 1:length(s)
   s(is).Color = "k"
   s(is).LineWidth = 1.3
   end
   =#
end

function dist_plot(pType='e')
   #dir = "/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia"
   dir = "/Users/hyzhou/"
   nBox = 4
   nbin = 60
   fs = 10
   subplotlist = [1,2,3,4,1,2,3,4]
   PlotVType = (1,1,1,1,3,3,3,3)

   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   fnameField = "3d_var_region0_0_"*fnameE[end-22:end]

   binRangeI = [[-3.,3.], [-3.,3.]]
   binRangeE = [[-7.,13.], [-10.,10.]]

   # Classify particles based on locations
   region = Array{Float32,2}(undef,6,nBox)
   particle = [Array{Float32}(undef, 3, 0) for _ in 1:nBox]

   # Electron
   if pType == 'e'
      region[:,1] = [-1.950, -1.945, -0.08, 0.08, -0.10, -0.06]
      region[:,2] = [-1.915, -1.910, -0.08, 0.08, -0.10, -0.06]
      region[:,3] = [-1.905, -1.900, -0.08, 0.08, -0.10, -0.06]
      region[:,4] = [-1.895, -1.890, -0.08, 0.08, -0.10, -0.06]

      head, data = readdata(fnameE, dir=dir)

      x = @view data[1].x[:,:,:,1]
      y = @view data[1].x[:,:,:,2]
      z = @view data[1].x[:,:,:,3]

      ux_ = findfirst(x->x=="ux", head[1][:wnames])
      uy_ = findfirst(x->x=="uy", head[1][:wnames])
      uz_ = findfirst(x->x=="uz", head[1][:wnames])

      ux = @view data[1].w[:,:,:,ux_]
      uy = @view data[1].w[:,:,:,uy_]
      uz = @view data[1].w[:,:,:,uz_]

      for ip = 1:length(x)
         for iR = 1:nBox
            if region[1,iR] < x[ip] < region[2,iR] &&
               region[3,iR] < y[ip] < region[4,iR] &&
               region[5,iR] < z[ip] < region[6,iR]

               particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
               break
            end
         end
      end

      binRange = binRangeE
   elseif pType == 'i'
      region[:,1] = [-1.865, -1.855, -0.08, 0.08, -0.29, -0.25]
      region[:,2] = [-1.910, -1.900, -0.08, 0.08, 0.11, 0.15]
      region[:,3] = [-1.970, -1.960, -0.08, 0.08, -0.10, -0.06]
      region[:,4] = [-1.850, -1.840, -0.08, 0.08, -0.10, -0.06]

      head, data = readdata(fnameI, dir=dir)

      x = @view data[1].x[:,:,:,1]
      y = @view data[1].x[:,:,:,2]
      z = @view data[1].x[:,:,:,3]

      ux_ = findfirst(x->x=="ux", head[1][:wnames])
      uy_ = findfirst(x->x=="uy", head[1][:wnames])
      uz_ = findfirst(x->x=="uz", head[1][:wnames])

      ux = @view data[1].w[:,:,:,ux_]
      uy = @view data[1].w[:,:,:,uy_]
      uz = @view data[1].w[:,:,:,uz_]

      for ip = 1:length(x)
         for iR = 1:nBox
            if region[1,iR] < x[ip] < region[2,iR] &&
               region[3,iR] < y[ip] < region[4,iR] &&
               region[5,iR] < z[ip] < region[6,iR]

               particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
               break
            end
         end
      end

      binRange = binRangeI
   end

   # Normalized quantities
   fig, ax = plt.subplots(4,2,figsize=(4.2,8.0))
   if pType == 'e'
      fig.suptitle("electron", fontsize=14, fontweight="bold")
   elseif pType == 'i'
      fig.suptitle("ion", fontsize=14, fontweight="bold")
   end

   xlPos = (0.50, -0.10)
   ylPos = (-0.25, 0.50)
   c = Vector{PyObject}(undef,nBox)
   plt.set_cmap("jet")

   for i = 1:nBox*2
      ax[i].tick_params(which="both", direction="in", top=true, right=true)
      ax[i].minorticks_on()
      iB = subplotlist[i]
      if PlotVType[i] ≤ 3
         ux = particle[iB][1,:] ./ cAlfven
         uy = particle[iB][2,:] ./ cAlfven
         uz = particle[iB][3,:] ./ cAlfven
      else
         dBx, dBy, dBz = GetMeanField(fnameField, region[:,iB]; dir=dir)

         dPar   = [dBx; dBy; dBz]         # Parallel direction
         dPerpI = cross([0; -1; 0], dPar) # Perpendicular direction in-plane
         dPerpO = cross(dPar, dPerpI)     # Perpendicular direction out-of-plane

         uPar   = transpose(particle[iB][1:3,:])*dPar   ./ cAlfven
         uPerpI = transpose(particle[iB][1:3,:])*dPerpI ./ cAlfven
         uPerpO = transpose(particle[iB][1:3,:])*dPerpO ./ cAlfven
      end

      if PlotVType[i]==1
         h = ax[i].hist2d(uy, ux, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[i]==2
         h = ax[i].hist2d(ux, uz, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[i]==3
         h = ax[i].hist2d(uy, uz, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[i]==4
         h = ax[i].hist2d(uPerpO, uPerpI, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[i]==5
         h = ax[i].hist2d(uPerpI, uPar, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[i]==6
         h = ax[i].hist2d(uPerpO, uPar, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      else
         @error "Unknown PlotVType!"
      end
      ax[i].set_aspect("equal", "box")
      #ax[i].grid(true)
      ax[i].axvline(x=0.0, ls="--", linewidth=0.5, color="k")
      ax[i].axhline(y=0.0, ls="--", linewidth=0.5, color="k")
      (i % 4 != 0) && ax[i].axes.xaxis.set_ticklabels([])

      if PlotVType[i]==1
         if i in (4,8)
            ax[i].annotate(L"u_y", xy=(0.50, -0.20),xycoords="axes fraction",
               fontsize=fs)
         else
            ax[i].annotate(L"u_y", xy=xlPos,xycoords="axes fraction",
               fontsize=fs)
         end
         ax[i].annotate(L"u_x", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[i]==2
         ax[i].annotate(L"u_x", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[i].annotate(L"u_z", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[i]==3
         if i in (4,8)
            ax[i].annotate(L"u_y", xy=(0.50, -0.20),xycoords="axes fraction",
               fontsize=fs)
         else
            ax[i].annotate(L"u_y", xy=xlPos,xycoords="axes fraction",
               fontsize=fs)
         end
         ax[i].annotate(L"u_z", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[i]==4
         ax[i].annotate(L"u_{\perp Out}", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[i].annotate(L"u_{\perp In}", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[i]==5
         ax[i].annotate(L"u_{\perp In}", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[i].annotate(L"u_\parallel", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[i]==6
         ax[i].annotate(L"u_{\perp Out}", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[i].annotate(L"u_\parallel", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      end

      if pType == 'e'
         ax[i].text(-0.05,1.02, @sprintf("x[%3.3f,%3.3f], z[%3.3f,%3.3f]",
            region[1,iB],region[2,iB],region[5,iB],region[6,iB]), FontSize=7,
            transform=ax[i].transAxes)
         ax[i].text(0.05,0.05,string(iB), FontSize=fs, color="r",
            transform=ax[i].transAxes)
      elseif pType == 'i'
         ax[i].text(-0.05,1.02, @sprintf("x[%3.3f,%3.3f], z[%3.3f,%3.3f]",
            region[1,iB],region[2,iB],region[5,iB],region[6,iB]),
            FontSize=7, transform=ax[i].transAxes)
         ax[i].text(0.83,0.05,string(iB+4), FontSize=fs, color="c",
            transform=ax[i].transAxes)
      end
   end

   im = plt.gca().get_children()[1]
   cax = fig.add_axes([0.1,0.93,0.8,0.02])
   colorbar(im, ax=ax, cax=cax, orientation="horizontal")
   cax.tick_params(which="both", axis="x", direction="out",color="k")

   return ax
end

function show_box_region()
   dir = "/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia"
   #dir = "/Users/hyzhou/"
   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   fnameField = "3d_var_region0_0_"*fnameE[end-22:end]

   fs = 10
   plotrange = [-2.0, -1.8, -0.35, 0.35]
   cutPlane = 129

   head, data = readdata(fnameField, dir=dir)

   X, Z, Bx = cutdata(data[1],head[1],"Bx",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)
   X, Z, Bz = cutdata(data[1],head[1],"Bz",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)
   X, Z, Ex = cutdata(data[1],head[1],"Ex",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)

   fig, ax = plt.subplots(1,1,figsize=(8.0,4.0))
   #plt.set_cmap("RdBu_r")
   plt.set_cmap("seismic")

   region = Array{Float32,2}(undef,6,8)
   region[:,1] = [-1.950, -1.945, -0.08, 0.08, -0.10, -0.06]
   region[:,2] = [-1.915, -1.910, -0.08, 0.08, -0.10, -0.06]
   region[:,3] = [-1.905, -1.900, -0.08, 0.08, -0.10, -0.06]
   region[:,4] = [-1.895, -1.890, -0.08, 0.08, -0.10, -0.06]

   region[:,5] = [-1.865, -1.855, -0.08, 0.08, -0.29, -0.25]
   region[:,6] = [-1.910, -1.900, -0.08, 0.08, 0.11, 0.15]
   region[:,7] = [-1.970, -1.960, -0.08, 0.08, -0.10, -0.06]
   region[:,8] = [-1.850, -1.840, -0.08, 0.08, -0.10, -0.06]

   c = ax.contourf(Z,X,Ex./1e3, 100, norm=matplotlib.colors.DivergingNorm(0),
      vmin=-12e1, vmax=12e1)


   # create an axes on the right side of ax. The width of cax will be 5%
   # of ax and the padding between cax and ax will be fixed at 0.05 inch.
   divider = axes_grid1.make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   colorbar(c, cax=cax)
   ax.invert_yaxis()
   ax.set_aspect("equal", "box")

   ax.set_xlabel(L"z [R_G]", fontsize=14)
   ax.set_ylabel(L"x [R_G]", fontsize=14)
   ax.set_title(L"Ex [mV/m]", fontsize=14)

   x  = @view X[:,1]
   z  = @view Z[1,:]

   seeds = select_seeds(x[10:end-10],z[10:end-10]; nSeed=1)
   xstart, zstart = seeds[1,:], seeds[2,:]
   append!(xstart, [-1.87, -1.97, -1.93, -1.93, -1.93, -1.95, -1.82, -1.93, -1.93, -1.93])
   append!(zstart, [-0.1, -0.1, -0.15, 0.0, 0.2, -0.25, -0.1, -0.08, -0.095, 0.116])

   xl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   zl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   for i = 1:length(xstart)
      xs, zs = xstart[i], zstart[i]
      xl[i], zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds=0.02, maxstep=20000,
         gridType="ndgrid")
   end
   [ax.plot(zl[j],xl[j],"-",color="k",lw=1.0) for j in 1:length(xstart)]

   ax.contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)

   for iB = 1:4
      rect = matplotlib.patches.Rectangle( (region[5,iB], region[1,iB]),
         region[6,iB]-region[5,iB], region[2,iB]-region[1,iB],
         lw=1.0, fill=true, alpha=0.9, fc="r") # facecolor="none"
      ax.add_patch(rect)
   end

   for iB = 5:8
      rect = matplotlib.patches.Rectangle( (region[5,iB], region[1,iB]),
         region[6,iB]-region[5,iB], region[2,iB]-region[1,iB],
         lw=1.0, fill=true, alpha=0.9, fc="c") # facecolor="none"
      ax.add_patch(rect)
   end

   ax.annotate("1", (region[5,1]-2e-2,region[2,1]), fontsize=10, color="r")
   ax.annotate("2", (region[5,2]-2e-2,region[2,2]), fontsize=10, color="r")
   ax.annotate("3", (region[5,3]-2e-2,region[2,3]+0.5e-2),fontsize=10,color="r")
   ax.annotate("4", (region[5,4]-2e-2,region[2,4]+1e-2), fontsize=10,color="r")

   for i = 5:8
      ax.annotate(string(i), (region[5,i]-2e-2, region[2,i]+0.2e-2), fontsize=10,
         color="c")
   end

   return ax
end

function HF_velocity()
   #dir = "/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia"
   dir = "/Users/hyzhou/"
   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   fnameField = "3d_var_region0_0_"*fnameE[end-22:end]

   head, data = readdata(fnameField, dir=dir)

   x = data[1].x[:,:,:,1]
   y = data[1].x[:,:,:,2]
   z = data[1].x[:,:,:,3]

   ux_ = findfirst(x->x=="uxs1", lowercase.(head[1][:wnames]))
   uy_ = findfirst(x->x=="uys1", lowercase.(head[1][:wnames]))
   uz_ = findfirst(x->x=="uzs1", lowercase.(head[1][:wnames]))

   Ux = @view data[1].w[:,:,:,ux_]
   Uy = @view data[1].w[:,:,:,uy_]
   Uz = @view data[1].w[:,:,:,uz_]

   #region[:,5] = [-1.850, -1.840, -0.08, 0.08, -0.29, -0.25]
   #region[:,6] = [-1.910, -1.900, -0.08, 0.08, 0.11, 0.15]
   limits = [-1.910, -1.900, -0.08, 0.08, 0.11, 0.15]
   xnew, ynew, znew, UxNew, UyNew, UzNew = subvolume(x,y,z, Ux,Uy,Uz, limits)

   # Average over the selected volume
   Ūx, Ūy, Ūz = mean(UxNew), mean(UyNew), mean(UzNew)

   return Ūx, Ūy, Ūz
end

#=
dir = "/Users/hyzhou"
fnameField = "3d_var_region0_0_t00001640_n00020369.out"
PType = 'e'
PlotVType = 2

if PType == 'e'
   fnameParticle = "cut_particles0_region0_1_t00001640_n00020369.out"
elseif PType == 'i'
   fnameParticle = "cut_particles1_region0_2_t00001640_n00020369.out"
end


# Define regions
xC, yC, zC = -1.75, 0.0, -0.2
xL, yL, zL = 0.008, 0.2, 0.03 # box length in x,y,z
=#

#=
@time region, particle = dist_select(
   fnameParticle, xC, yC, zC, xL, yL, zL,
   dir=dir, ParticleType=PType)

@time dist_scan(region, particle, PType, PlotVType; dir=dir, fnameField=fnameField)

@time plotExCut(fnameField, region, xC,yC,zC,xL,yL,zL, dir=dir)
=#

#ax = dist_plot('i')
#ux, uy, uz = HF_velocity()
show_box_region() # horizontal
