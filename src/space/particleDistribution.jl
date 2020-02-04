# Distribution plot in a large PIC domain.
#
#
# Hongyang Zhou, hyzhou@umich.edu 01/30/2020

using VisAna, PyPlot, Printf

## Parameters
const cAlfven = 253       # average Alfven speed in G8, [km/s]
const me = 9.10938356e-31 # electron mass, [kg]
const mp = 1.6726219e-27  # proton mass, [kg]
const mi = 14             # average ion mass [amu]
const nBox = 9            # number of box regions

"""
	dist_select(fnameParticle; ParticleType='e', dir=".")

Select particle in regions.
`ParticleType` in ['e','i'].
xC, yC, zC are the central box center position.
"""
function dist_select(fnameParticle, xC=-1.90, yC=0.0, zC=-0.1; ParticleType='e',
   dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia/",)

   if ParticleType == 'i'
      !occursin("region0_2", fnameParticle) && @error "Check filename!"
   elseif ParticleType == 'e'
      !occursin("region0_1", fnameParticle) && @error "Check filename!"
   end

   fnameField = "3d_var_region0_0_"*fnameParticle[end-22:end]

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
PlotVType: 1: v_par vs. v_perp1 2: v_perp1 vs. v_perp2
"""
function dist_plot(region, particle, ParticleType='i', PlotVType=1)

   if ParticleType == 'i'
      binRange = [[-3.,3.], [-3.,3.]]
   elseif ParticleType == 'e'
      binRange = [[-10.,12.], [-10.,10.]]
   end

   figure(figsize=(11,6))
   for iB = 1:nBox
      ux = particle[iB][1,:] ./ cAlfven
      uy = particle[iB][2,:] ./ cAlfven
      uz = particle[iB][3,:] ./ cAlfven

      ax = subplot(3,4,iB+ceil(iB/3))
      if PlotVType==1
         h = hist2D(uy, ux, bins=60,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==2
         h = hist2D(ux, uz, bins=60,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      elseif PlotVType==3
         h = hist2D(uy, uz, bins=60,
            norm=matplotlib.colors.LogNorm(),density=true, range=binRange)
      else
         @error "Unknown PlotVType!"
      end
      grid(true)
      axis("equal")

      if PlotVType==1
         xlabel(L"u_y",fontsize=14)
         ylabel(L"u_x",fontsize=14)
      elseif PlotVType==2
         xlabel(L"u_x",FontSize=14)
         ylabel(L"u_z",FontSize=14)
      elseif PlotVType==3
         xlabel(L"u_y",FontSize=14)
         ylabel(L"u_x",FontSize=14)
      end
      title(@sprintf("%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]",iB,region[1,iB],
         region[2,iB],region[5,iB],region[6,iB]))
      colorbar()
      plt.set_cmap("hot")

      if ParticleType == 'e'
         str = "electron"
      elseif ParticleType == 'i'
         str = "ion"
      end
      text(0.05,0.05,str, FontSize=14, transform=ax.transAxes)
   end

end



function plotExCut(fnameField::String, region, xC, yC, zC, xL, yL, zL;
   dir="/Users/hyzhou")

   plotrange = [xC-xL*16, xC+xL*16, zC-zL*5, zC+zL*5]
   # Sample region plot over contour
   @time head, data = readdata(fnameField, dir=dir)

   bx_ = findfirst(x->x=="Bx", head[1][:wnames])
   bz_ = findfirst(x->x=="Bz", head[1][:wnames])
   #ex_ = findfirst(x->x=="Ex", head[1][:wnames])

   Bx = @view data[1].w[:,:,:,bx_]
   Bz = @view data[1].w[:,:,:,bz_]
   #Ex = @view data[1].w[:,:,:,ex_]

   subplot(3,4,(1,9))
   cutplot(data[1],head[1],"Ex",cut='y',cutPlaneIndex=128,plotrange=plotrange)
   colorbar()
   axis("scaled")
   plt.set_cmap("RdBu_r")
   clim(-9e4,9e4)
   xlabel(L"x [R_G]", fontsize=16)
   ylabel(L"z [R_G]", fontsize=16)
   title(L"Ex [\mu V/m]")

   for iB = 1:nBox
      rect = matplotlib.patches.Rectangle( (region[1,iB], region[5,iB]),
      region[2,iB]-region[1,iB], region[6,iB]-region[5,iB],
      ec="r", lw=1.2, fill=false) # facecolor="none"
      ax = gca()
      #ax = subplot(3,4,(1,9))
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
tight_layout()

end


PType = 'e'
PlotVType = 1
# Define regions
xC, yC, zC = -1.90, 0.0, -0.1
xL, yL, zL = 0.005, 0.2, 0.07 # box length in x,y,z

@time region, particle = dist_select(
   "cut_particles0_region0_1_t00001640_n00020369.out", xC, yC, zC, xL, yL, zL,
   dir="/Users/hyzhou", ParticleType=PType)
@time dist_plot(region, particle, PType, PlotVType)

@time plotExCut("3d_var_region0_0_t00001640_n00020369.out", region,
   xC,yC,zC,xL,yL,zL, dir="/Users/hyzhou")
