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
	dist_plot(region, particle, ParticleType='i', PlotVType=1)

Velocity distribution plot in 9 regions.
`PlotVType`: 1: uy-ux; 2: ux-uz; 3:uy-uz; 4:u⟂O-u⟂I; 5:u⟂I-u∥; 5:u⟂O-u∥
"""
function dist_plot(region, particle, ParticleType='i', PlotVType=1; dir=".",
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

function NicePlot()
   dir = "/Users/hyzhou"
   fnameField = "3d_var_region0_0_t00001640_n00020369.out"
   PlotVType = 1
   nBox = 8
   nbin = 60
   fs = 10

   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   fnameField = "3d_var_region0_0_"*fnameE[end-22:end]

   # Classify particles based on locations
   region = Array{Float32,2}(undef,6,nBox)
   region[:,1] = [-1.930, -1.925, -0.08, 0.08, -0.10, -0.06]
   region[:,2] = [-1.915, -1.910, -0.08, 0.08, -0.10, -0.06]
   region[:,3] = [-1.905, -1.900, -0.08, 0.08, -0.10, -0.06]
   region[:,4] = [-1.890, -1.885, -0.08, 0.08, -0.10, -0.06]

   region[:,5] = [-1.905, -1.900, -0.08, 0.08, -0.10, -0.06]
   region[:,6] = [-1.900, -1.895, -0.08, 0.08, -0.10, -0.06]
   region[:,7] = [-1.895, -1.890, -0.08, 0.08, -0.10, -0.06]
   region[:,8] = [-1.890, -1.885, -0.08, 0.08, -0.10, -0.06]

   particle = [Array{Float32}(undef, 3, 0) for _ in 1:nBox]

   # Electron
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
      for iR = 1:Int(nBox/2)
         if region[1,iR] < x[ip] < region[2,iR] &&
            region[3,iR] < y[ip] < region[4,iR] &&
            region[5,iR] < z[ip] < region[6,iR]

            particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
            break
         end
      end
   end

   # Ion
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
      for iR = Int(nBox/2)+1:nBox
         if region[1,iR] < x[ip] < region[2,iR] &&
            region[3,iR] < y[ip] < region[4,iR] &&
            region[5,iR] < z[ip] < region[6,iR]

            particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
            break
         end
      end
   end


   binRangeI = [[-3.,3.], [-3.,3.]]
   binRangeE = [[-10.,12.], [-10.,10.]]

   figure(figsize=(4,10))
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

      ax = subplot(5,2,2+iB)
      if iB ≤ nBox/2
         binRange = binRangeE
      else
         binRange = binRangeI
      end
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
      ax.set_aspect("equal", "box")
      #axis("equal")

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
      title(@sprintf("%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]",iB,region[1,iB],
         region[2,iB],region[5,iB],region[6,iB]))
      colorbar()
      plt.set_cmap("jet")
      #clim(1e-2,10^0.3)

      if iB ≤ 4
         text(0.05,0.05,"electron", FontSize=fs, transform=ax.transAxes)
      else
         text(0.05,0.05,"ion", FontSize=fs, transform=ax.transAxes)
      end
   end

   plotrange = [-2.05, -1.75, -0.5, 0.5]
   cutPlane = 129

   head, data = readdata(fnameField, dir=dir)

   X, Z, Bx = cutdata(data[1],head[1],"Bx",cut='y',cutPlaneIndex=cutPlane,
   	plotrange=plotrange)
   X, Z, Bz = cutdata(data[1],head[1],"Bz",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)

   ax = subplot(5,2,(1,2))
   X, Z, Ex = cutdata(data[1],head[1],"Ex",cut='y',cutPlaneIndex=cutPlane,
   	plotrange=plotrange)
   c = contourf(Z,X,Ex, 40, norm=matplotlib.colors.DivergingNorm(0), vmin=-12e4,
      vmax=12e4)

   # create an axes on the right side of ax. The width of cax will be 5%
   # of ax and the padding between cax and ax will be fixed at 0.05 inch.
   divider = axes_grid1.make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   colorbar(c, cax=cax)
   ax.invert_yaxis()
   ax.set_aspect("equal", "box")
   plt.set_cmap("RdBu_r")
   #clim(-9e4,9e4)

   ax.set_xlabel(L"z [R_G]", fontsize=fs)
   ax.set_ylabel(L"x [R_G]", fontsize=fs)
   ax.set_title(L"Ex [\mu V/m]")

   x  = @view X[:,1]
   z  = @view Z[1,:]
   zstart = collect(range(z[10],stop=z[end-10],length=14))
   xstart = fill(-1.92,size(zstart))
   zstart = append!(zstart, collect(range(-0.3,0.4,length=3)))
   xstart = append!(xstart, fill(-1.82,3))

   xl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   zl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   for i = 1:length(xstart)
      xs,zs = xstart[i],zstart[i]
      xl[i], zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds=0.02, maxstep=20000,
      gridType="ndgrid")
   end
   [ax.plot(zl[j],xl[j],"-",color="k",lw=0.4) for j in 1:length(xstart)]

   for iB = 1:nBox
      if iB ≤ nBox/2
         rect = matplotlib.patches.Rectangle( (region[5,iB], region[1,iB]),
            region[6,iB]-region[5,iB], region[2,iB]-region[1,iB],
            ec="r", lw=1.2, fill=false) # facecolor="none"
      else
         rect = matplotlib.patches.Rectangle( (region[5,iB], region[1,iB]),
            region[6,iB]-region[5,iB], region[2,iB]-region[1,iB],
            ec="b", lw=1.2, fill=false) # facecolor="none"
      end
      ax.add_patch(rect)
   end

end

function NicePlot2()
   dir = "/Users/hyzhou"
   nBox = 8
   nbin = 60
   fs = 10

   fnameE = "cut_particles0_region0_1_t00001640_n00020369.out"
   fnameI = "cut_particles1_region0_2_t00001640_n00020369.out"

   fnameField = "3d_var_region0_0_"*fnameE[end-22:end]

   # Classify particles based on locations
   region = Array{Float32,2}(undef,6,nBox)
   region[:,1] = [-1.930, -1.925, -0.08, 0.08, -0.10, -0.06]
   region[:,2] = [-1.915, -1.910, -0.08, 0.08, -0.10, -0.06]
   region[:,3] = [-1.905, -1.900, -0.08, 0.08, -0.10, -0.06]
   region[:,4] = [-1.890, -1.885, -0.08, 0.08, -0.10, -0.06]

   region[:,5] = [-1.905, -1.900, -0.08, 0.08, -0.10, -0.06]
   region[:,6] = [-1.900, -1.895, -0.08, 0.08, -0.10, -0.06]
   region[:,7] = [-1.895, -1.890, -0.08, 0.08, -0.10, -0.06]
   region[:,8] = [-1.890, -1.885, -0.08, 0.08, -0.10, -0.06]

   particle = [Array{Float32}(undef, 3, 0) for _ in 1:nBox]

   # Electron
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
      for iR = 1:Int(nBox/2)
         if region[1,iR] < x[ip] < region[2,iR] &&
            region[3,iR] < y[ip] < region[4,iR] &&
            region[5,iR] < z[ip] < region[6,iR]

            particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
            break
         end
      end
   end

   # Ion
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
      for iR = Int(nBox/2)+1:nBox
         if region[1,iR] < x[ip] < region[2,iR] &&
            region[3,iR] < y[ip] < region[4,iR] &&
            region[5,iR] < z[ip] < region[6,iR]

            particle[iR] = hcat(particle[iR], [ux[ip]; uy[ip]; uz[ip]])
            break
         end
      end
   end


   binRangeI = [[-3.,3.], [-3.,3.]]
   binRangeE = [[-10.,12.], [-10.,10.]]

   # Normalized quantities
   fig, ax = plt.subplots(4,3,figsize=(10.0,10.0))
   xlPos = (0.5, -0.12)
   ylPos = (-0.15, 0.5)
   subplotlist = (5,6,9,10,7,8,11,12)
   PlotVType = (1,1,3,3,1,1,3,3)
   c = Vector{PyObject}(undef,nBox)
   axin = Vector{PyObject}(undef,nBox)
   for iB in 1:nBox
   	axin[iB] = inset_axes(ax[subplotlist[iB]],
         width="5%",  # width = 5% of parent_bbox width
         height="100%",  # height : 50%
         loc="lower left",
         bbox_to_anchor=(1.02, 0., 1.0, 1.0),
         bbox_transform=ax[subplotlist[iB]].transAxes,
         borderpad=0,)
      axin[iB].tick_params(axis="y", direction="in")
      
      ax[subplotlist[iB]].tick_params(which="both", direction="in",
         top=true, right=true)
   end


   plt.set_cmap("jet")
   for iB = 1:nBox
      if PlotVType[iB] ≤ 3
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

      #ax = subplot(4,3,subplotlist[iB])

      if iB ≤ nBox/2
         binRange = binRangeE
      else
         binRange = binRangeI
      end
      if PlotVType[iB]==1
         h = ax[subplotlist[iB]].hist2d(uy, ux, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[iB]==2
         h = ax[subplotlist[iB]].hist2d(ux, uz, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[iB]==3
         h = ax[subplotlist[iB]].hist2d(uy, uz, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[iB]==4
         h = ax[subplotlist[iB]].hist2d(uPerpO, uPerpI, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[iB]==5
         h = ax[subplotlist[iB]].hist2d(uPerpI, uPar, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType[iB]==6
         h = ax[subplotlist[iB]].hist2d(uPerpO, uPar, bins=nbin,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      else
         @error "Unknown PlotVType!"
      end
      colorbar(h[4], cax=axin[iB])
      #clim(1e-2,10^0.3)
      ax[subplotlist[iB]].set_aspect("equal", "box")
      #axis("equal")
      ax[subplotlist[iB]].grid(true)

      if PlotVType[iB]==1
         ax[subplotlist[iB]].annotate(L"u_y", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[subplotlist[iB]].annotate(L"u_x", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[iB]==2
         ax[subplotlist[iB]].annotate(L"u_x", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[subplotlist[iB]].annotate(L"u_z", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[iB]==3
         ax[subplotlist[iB]].annotate(L"u_y", xy=xlPos,xycoords="axes fraction",
            fontsize=fs)
         ax[subplotlist[iB]].annotate(L"u_z", xy=ylPos,xycoords="axes fraction",
            fontsize=fs)
      elseif PlotVType[iB]==4
         ax[subplotlist[iB]].set_xlabel(L"u_{\perp Out}",fontsize=fs)
         ax[subplotlist[iB]].set_ylabel(L"u_{\perp In}",fontsize=fs)
      elseif PlotVType[iB]==5
         ax[subplotlist[iB]].set_xlabel(L"u_{\perp In}",FontSize=fs)
         ax[subplotlist[iB]].set_ylabel(L"u_\parallel",FontSize=fs)
      elseif PlotVType[iB]==6
         ax[subplotlist[iB]].set_xlabel(L"u_{\perp Out}",FontSize=fs)
         ax[subplotlist[iB]].set_ylabel(L"u_\parallel",FontSize=fs)
      end
      ax[subplotlist[iB]].set_title(@sprintf("%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]",iB,region[1,iB],
         region[2,iB],region[5,iB],region[6,iB]), FontSize=5)

      if iB ≤ 4
         ax[subplotlist[iB]].text(0.05,0.05,"electron", FontSize=fs,
            transform=ax[subplotlist[iB]].transAxes)
      else
         ax[subplotlist[iB]].text(0.05,0.05,"ion", FontSize=fs,
            transform=ax[subplotlist[iB]].transAxes)
      end
   end

   plotrange = [-2.05, -1.75, -0.5, 0.5]
   cutPlane = 129

   head, data = readdata(fnameField, dir=dir)

   X, Z, Bx = cutdata(data[1],head[1],"Bx",cut='y',cutPlaneIndex=cutPlane,
   	plotrange=plotrange)
   X, Z, Bz = cutdata(data[1],head[1],"Bz",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)
   X, Z, Ex = cutdata(data[1],head[1],"Ex",cut='y',cutPlaneIndex=cutPlane,
   	plotrange=plotrange)

   plt.set_cmap("RdBu_r")
   ax = subplot(4,3,(1,4))
   c = contourf(X,Z,Ex./1e3, 40, norm=matplotlib.colors.DivergingNorm(0), vmin=-12e1,
      vmax=12e1)

   # create an axes on the right side of ax. The width of cax will be 5%
   # of ax and the padding between cax and ax will be fixed at 0.05 inch.
   divider = axes_grid1.make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   colorbar(c, cax=cax)
   #ax.invert_yaxis()
   ax.set_aspect("equal", "box")
   #clim(-9e4,9e4)

   ax.set_xlabel(L"z [R_G]", fontsize=fs)
   ax.set_ylabel(L"x [R_G]", fontsize=fs)
   ax.set_title(L"Ex [mV/m]")

   x  = @view X[:,1]
   z  = @view Z[1,:]
   zstart = collect(range(z[10],stop=z[end-10],length=14))
   xstart = fill(-1.92,size(zstart))
   zstart = append!(zstart, collect(range(-0.3,0.4,length=3)))
   xstart = append!(xstart, fill(-1.82,3))

   xl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   zl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   for i = 1:length(xstart)
      xs,zs = xstart[i],zstart[i]
      xl[i], zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds=0.02, maxstep=20000,
      gridType="ndgrid")
   end
   [ax.plot(xl[j],zl[j],"-",color="k",lw=0.4) for j in 1:length(xstart)]

   for iB = 1:Int(nBox/2)
      rect = matplotlib.patches.Rectangle( (region[1,iB], region[5,iB]),
         region[2,iB]-region[1,iB], region[6,iB]-region[5,iB],
         ec="r", lw=1.2, fill=false) # facecolor="none"
      ax.add_patch(rect)
   end

   ax = subplot(4,3,(7,10))
   X, Z, Ex = cutdata(data[1],head[1],"Ex",cut='y',cutPlaneIndex=cutPlane,
      plotrange=plotrange)
   c = contourf(X,Z,Ex./1e3, 40, norm=matplotlib.colors.DivergingNorm(0), vmin=-12e1,
      vmax=12e1)

   # create an axes on the right side of ax. The width of cax will be 5%
   # of ax and the padding between cax and ax will be fixed at 0.05 inch.
   divider = axes_grid1.make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   colorbar(c, cax=cax)
   #ax.invert_yaxis()
   ax.set_aspect("equal", "box")
   plt.set_cmap("RdBu_r")
   #clim(-9e4,9e4)

   ax.set_xlabel(L"z [R_G]", fontsize=fs)
   ax.set_ylabel(L"x [R_G]", fontsize=fs)
   ax.set_title(L"Ex [mV/m]")

   x  = @view X[:,1]
   z  = @view Z[1,:]
   zstart = collect(range(z[10],stop=z[end-10],length=14))
   xstart = fill(-1.92,size(zstart))
   zstart = append!(zstart, collect(range(-0.3,0.4,length=3)))
   xstart = append!(xstart, fill(-1.82,3))

   xl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   zl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
   for i = 1:length(xstart)
      xs,zs = xstart[i],zstart[i]
      xl[i], zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds=0.02, maxstep=20000,
      gridType="ndgrid")
   end
   [ax.plot(xl[j],zl[j],"-",color="k",lw=0.4) for j in 1:length(xstart)]

   for iB = Int(nBox/2)+1:nBox
      rect = matplotlib.patches.Rectangle( (region[1,iB], region[5,iB]),
         region[2,iB]-region[1,iB], region[6,iB]-region[5,iB],
         ec="b", lw=1.2, fill=false) # facecolor="none"
      ax.add_patch(rect)
   end

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

@time dist_plot(region, particle, PType, PlotVType; dir=dir, fnameField=fnameField)

@time plotExCut(fnameField, region, xC,yC,zC,xL,yL,zL, dir=dir)
=#

NicePlot2()
