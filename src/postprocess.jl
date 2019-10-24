# Script for further analysis of simulation results.
#
# Hongyang Zhou, hyzhou@umich.edu 10/01/2019

using Pkg
Pkg.activate("/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia")

using VisAna, PyPlot, PyCall, Glob, Printf

include("constants.jl")

"""Series of contour plots from box output."""
function plot_contour_from_box(filename::String, dir::String=".")

   filehead, data, filelist = readdata(filename, dir=dir, verbose=false);
   npict = filelist[1].npictinfiles

   saveDir = "ycutsFromBox"
   mkdir(saveDir)

   for ipict = 1:npict
      filehead, data, filelist = readdata(filename, npict=ipict, dir=dir,
         verbose=false)

      plotdata(data[1],filehead[1],"by bx;bz",plotmode="contbar streamover",
         cut="y",cutPlaneIndex=3, density=1.5)

      fig = matplotlib.pyplot.gcf()
      fig.set_size_inches(3.5,6)
      plt.axis("scaled")

      plt.savefig(saveDir*"/y_"*string(ipict)*".png")
      close()
   end

end

"""Series of contour plots from cuts."""
function plot_contour_from_cuts(filename::String, dir::String=".")

   filehead, data, filelist = readdata(filename, dir=dir, verbose=false);
   npict = filelist[1].npictinfiles
   saveDir = "ycuts"
   mkdir(saveDir)

   for ipict = 1:npict
      filehead, data, filelist = readdata(filename,npict=ipict,dir=dir,
         verbose=false)

      plotdata(data[1], filehead[1], "by bx;bz", plotmode="contbar streamover",
         plotrange=[1.0 1.8 -2.0 2.0], density=1.5)

      fig = matplotlib.pyplot.gcf()
      fig.set_size_inches(3.5,6)
      plt.axis("scaled")

      plt.savefig(saveDir*"/y_"*string(ipict)*".png")
      close()
   end

end



"""
   get_diamagnetic_velocity(filename, filedir, nType)
Calculate diamagnetic velocity from outputs.

nType = 1 is Hall MHD, nType = 2 is PIC.
"""
function get_diamagnetic_velocity(filename::String, filedir::String=".",
   nType::Int=1)

   saveDir = "ycuts"
   mkdir(saveDir)

   if filename[end-3:end] == "outs"
      filehead, data, filelist = readdata(filename, dir=dir, verbose=false)
      npict = filelist[1].npictinfiles
      for i = 1:npict
         filehead, data, filelist =
            readdata(filename, dir=filedir, npict=i, verbose=false)
         processing(filehead[1], data[1], saveDir, nType)
      end
   else
      filenames = glob(filename, filedir)
      npict = length(filenames)
      for i = 1:npict
         filehead, data, filelist =
            readdata(filenames[i], dir=filedir, verbose=false)
         processing(filehead[1], data[1], saveDir, nType)
      end
   end

end

function processing(filehead::Dict, data::Data, saveDir::String, nType::Int=1)

   np = pyimport("numpy")

   @views x = data.x[:,:,:,1]
   @views y = data.x[:,:,:,2]
   @views z = data.x[:,:,:,3]
   w = data.w

   if nType == 1 # Hall MHD
      ni_ = findfirst(isequal("Rho"), filehead[:wnames])
      pi_ = findfirst(isequal("P"),   filehead[:wnames])
   elseif nType == 2 # PIC
      ni_ = findfirst(isequal("rhos1"), filehead[:wnames])
      pi_ = findfirst(isequal("ps1"),   filehead[:wnames])
   end

   bx_ = findfirst(isequal("Bx"),  filehead[:wnames])
   by_ = findfirst(isequal("By"),  filehead[:wnames])
   bz_ = findfirst(isequal("Bz"),  filehead[:wnames])

   ni = w[:,:,:,ni_] # [/cm^3]
   pi = w[:,:,:,pi_] # [nPa]
   Bx = w[:,:,:,bx_] # [nT]
   By = w[:,:,:,by_] # [nT]
   Bz = w[:,:,:,bz_] # [nT]

   # This assumes uniform grid resolution!
   px, py, pz = np.gradient(pi, x[2]-x[1])

   nX, nY, nZ = size(w)

   vDiaMag = Array{Float64,4}(undef,3,nX,nY,nZ)
   vMag    = Array{Float64,3}(undef,nX,nY,nZ)

   cUnit = 1.0e7/(14*1.602*2.634)

   # v = B×▽p/neB^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag2 = Bx[i,j,k]^2 + Bx[i,j,k]^2 + Bz[i,j,k]^2
      vDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) /
         (ni[i,j,k]*Bmag2) * cUnit
      vDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) /
         (ni[i,j,k]*Bmag2) * cUnit
      vDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) /
         (ni[i,j,k]*Bmag2) * cUnit
      vMag[i,j,k] = sqrt(
         vDiaMag[1,i,j,k]^2 + vDiaMag[2,i,j,k]^2 + vDiaMag[3,i,j,k]^2)
   end

   yMid = floor(Int,nY/2)
   xrange = floor(Int,nX/4):floor(Int,nX/4*3)
   zrange = floor(Int,nZ/4):floor(Int,nZ/4*3)

   clf()

   #cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:],levels=np.logspace(6,8,30))
   #cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:],
   #   locator=matplotlib.ticker.LogLocator(subs=collect(1:10)))
   #cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:], levels=50)
   #cf = contourf(x[xrange,yMid,zrange],z[xrange,yMid,zrange],vMag[xrange,yMid,zrange], levels=50)
   #cf = contourf(x[xrange,yMid,zrange],z[xrange,yMid,zrange],vMag[xrange,yMid,zrange], locator=matplotlib.ticker.LogLocator(), extend="both")
   v = np.linspace(0.0, 20.0, 11, endpoint=true)
   cf = contourf(x[xrange,yMid,zrange],z[xrange,yMid,zrange],log.(vMag[xrange,yMid,zrange]), v, levels=50, extend="both")

   fig = matplotlib.pyplot.gcf()
   fig.set_size_inches(3.5,6)
   plt.axis("scaled")
   cbar = plt.colorbar(cf,boundaries=v, ticks=v)
   plt.clim(0.0, 20.0)
   #plt.clim(0, 3e6) # This creates the same color for values above
   xlabel("x"); ylabel("z")
   title(L"$\mathbf{B}\times \nabla P/(neB^2)$"*", t=$(i)")
   plt.savefig(saveDir*"/diamagnetic_$(i).png")
end

function vtk_auto_conversion(fnames::String, dir::String=".")

   filenames = Vector{String}(undef,0)
   filesfound = glob(fnames, dir)
   filenames = vcat(filenames, filesfound)
   for filename in filenames
      head, data, connectivity = readtecdata(filename, false)
      convertVTK(head, data, connectivity, filename[1:end-4])
   end

end

#=
filename = "box*outs";
dir = "/Users/hyzhou/Ganymede/run_mercury_80s_newbox/GM";

filename = "y*outs";
dir = "/Users/hyzhou/Ganymede/run_mercury_80s_newbox/GM";
=#

#filename = "cut*.dat";
#dir = "GM/IO2";
#vtk_auto_conversion(filename, dir)


using Pkg
Pkg.activate("/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia");
using VisAna, PyCall, PyPlot

function y_by_cuts()
   np = pyimport("numpy");
   filename = "y*out";
   filehead, data, filelist = readdata(filename,verbose=false);

   fig, ax = plt.subplots(1, 1);
   var = "by"; plotinterval = 0.1;
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(filehead[1][:wnames]));
   X = vec(data[1].x[:,:,1]);
   Y = vec(data[1].x[:,:,2]);
   W = vec(data[1].w[:,:,VarIndex_]);

   # Create grid values first.
   xi = range(minimum(X), stop=maximum(X), step=plotinterval);
   yi = range(minimum(Y), stop=maximum(Y), step=plotinterval);
   # Perform linear interpolation of the data (x,y) on grid(xi,yi)
   triang = matplotlib.tri.Triangulation(X,Y);
   interpolator = matplotlib.tri.LinearTriInterpolator(triang, W);
   Xi, Yi = np.meshgrid(xi, yi);
   wi = interpolator(Xi, Yi);

   c = contourf(xi,yi,wi,50);
   fig.set_size_inches(6,6);
   plt.axis("scaled");
   colorbar(c, ax=ax,ticks=np.arange(-100.0, 0.0, 12.5));
   xlabel("x"); ylabel("z")
   title(L"By")
   plt.savefig("test.png")
end

function plot_beta(filename::String)

   filehead, data, filelist = readdata(filename,verbose=false);

   cutPlaneIndex = 65
   VarIndex_ = 18

   filehead = filehead[1]
   data = data[1]

   X = @view data.x[:,:,:,1]
   Y = @view data.x[:,:,:,2]
   Z = @view data.x[:,:,:,3]

   fig, ax = plt.subplots(1, 1)
   fig.set_size_inches(3.3,6);

   W = @view data.w[:,:,:,VarIndex_]
   Bx = @view data.w[:,:,:,5]
   By = @view data.w[:,:,:,6]
   Bz = @view data.w[:,:,:,7]

   cut1 = @view X[:,cutPlaneIndex,:]
   cut2 = @view Z[:,cutPlaneIndex,:]
   W    = @view W[:,cutPlaneIndex,:]
   Bx   = @view Bx[:,cutPlaneIndex,:]
   By   = @view By[:,cutPlaneIndex,:]
   Bz   = @view Bz[:,cutPlaneIndex,:]
   PB   = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)

   c = ax.contourf(cut1,cut2,W./PB)
   fig.colorbar(c,ax=ax)
   ax.axis("scaled")
   title(filehead[:wnames][VarIndex_])

   xlabel("x"); ylabel("z")

   dim = [0.125, 0.013, 0.2, 0.045]
   str = @sprintf "it=%d, time=%4.2f" filehead[:it] filehead[:time]
   at = matplotlib.offsetbox.AnchoredText(str,
      loc="lower left", prop=Dict("size"=>8), frameon=true,
      bbox_to_anchor=(0., 1.),
      bbox_transform=ax.transAxes)
      at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      ax.add_artist(at)
   #plt.savefig("beta.png")
end
