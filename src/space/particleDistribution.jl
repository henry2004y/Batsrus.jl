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
ParticleType: {'e','i'}
"""
function dist_select(fnameParticle; ParticleType='e',
   dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia/",)

   if ParticleType == 'i'
      fnameParticle = "cut_particles1_region0_2_t00001520_n00004093.out"
   elseif ParticleType == 'e'
      fnameParticle = "cut_particles0_region0_1_t00001520_n00004093.out"
   end

   fnameField = "3d_var_region0_0_"*fnameParticle[end-22:end]

   # Define regions
   xC = -1.80   # center of boxes
   yC = 0.0
   zC = -0.25
   xL = 0.005 # box length in x
   yL = 0.2   # box length in y
   zL = 0.07  # box length in z

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

   figure(figsize=(11,7))
   for ipict = 1:nBox
      ux = particle[ipict][1,:] ./ cAlfven
      uy = particle[ipict][2,:] ./ cAlfven
      uz = particle[ipict][3,:] ./ cAlfven

      ax = subplot(3,4,ipict+ceil(ipict/3))
      if PlotVType==1
         h = hist2D(uy, ux, bins=60,
            norm=matplotlib.colors.LogNorm(),range=[[-10.,12.], [-10.,12.]])
      elseif PlotVType==2
         h = hist2D(ux, uz, bins=60,
            norm=matplotlib.colors.LogNorm(),range=[[-3.,3.], [-3.,3.]])
      else
         @error "Unknown PlotVType!"
      end
      axis("equal")

      if PlotVType==1
         xlabel(L"u_y",fontsize=14)
         ylabel(L"u_x",fontsize=14)
      elseif PlotVType==2
         xlabel(L"u_x",FontSize=14)
         ylabel(L"u_z",FontSize=14)
      end
      title(@sprintf("%d, x[%3.3f,%3.3f], z[%3.3f,%3.3f]",ipict,region[1,ipict],
         region[1,ipict],region[5,ipict],region[6,ipict]))
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

ParticleType = 'e'
PlotVType = 1

@time region, particle = dist_select("cut_particles0_region0_1_t00001520_n00004093.out")
@time dist_plot(region, particle, ParticleType, PlotVType)

filename = "3d_var_region0_0_t00001520_n00004093.out"
# Sample region plot over contour
head, data = readdata(filename, dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia/")

# Define regions
xC = -1.80   # center of boxes
yC = 0.0
zC = -0.25
xL = 0.005 # box length in x
yL = 0.2   # box length in y
zL = 0.07  # box length in z
plotrange = [xC-xL*16, xC+xL*16, zC-zL*5, zC+zL*5]

bx_ = findfirst(x->x=="Bx", head[1][:wnames])
bz_ = findfirst(x->x=="Bz", head[1][:wnames])
ex_ = findfirst(x->x=="Ex", head[1][:wnames])

Bx = @view data[1].w[:,:,:,bx_]
Bz = @view data[1].w[:,:,:,bz_]
Ex = @view data[1].w[:,:,:,ex_]


subplot(3,4,(1,9))
cutplot(data[1],head[1],"Ex",cut='y',cutPlaneIndex=128,plotrange=plotrange)
colorbar()
axis("equal")
set_cmap("*RdBu")
#caxis([-9e4,9e4])
xlabel("x [R_G]", fontsize=16)
ylabel("z [R_G]", fontsize=16)
title(L"Ex [\mu V/m]")
#=
# streamline function requires the meshgrid format strictly
s = streamslice(cut1",cut2",Bx",Bz",1,"linear")
for is = 1:length(s)
   s(is).Color = "k"
   s(is).LineWidth = 1.3
end

for ipict = 1:nBox
   rectangle("Position",[Region{ipict}(1) Region{ipict}(5) ...
      Region{ipict}(2)-Region{ipict}(1) ...
      Region{ipict}(6)-Region{ipict}(5)],"EdgeColor","r","LineWidth",1.5)
end
=#
