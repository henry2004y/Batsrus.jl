# Script for further analysis of simulation results.
#
# Hongyang Zhou, hyzhou@umich.edu 10/01/2019

using VisAna, PyPlot, PyCall, Glob

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

   ni = w[:,:,:,ni_]*1e6  # [/m^3]
   pi = w[:,:,:,pi_]*1e-9 # [Pa]
   Bx = w[:,:,:,bx_]*1e-9 # [T]
   By = w[:,:,:,by_]*1e-9 # [T]
   Bz = w[:,:,:,bz_]*1e-9 # [T]

   # This assumes uniform grid resolution!
   px, py, pz = np.gradient(pi,(x[2]-x[1])*Rg)

   nX, nY, nZ = size(w)

   vDiaMag = Array{Float64,4}(undef,3,nX,nY,nZ)
   vMag    = Array{Float64,3}(undef,nX,nY,nZ)

   # v = B×▽p/neB^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag2 = Bx[i,j,k]^2 + Bx[i,j,k]^2 + Bz[i,j,k]^2
      vDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
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

function vtk_auto_conversion(fnames::String)

   dir = "."
   filenames = Vector{String}(undef,0)
   filesfound = glob(fnames, dir)
   filenames = vcat(filenames, filesfound)
   for filename in filenames
      head, data, connectivity = readtecdata(filename, false)
      convertVTK(head, data, connectivity, filename[1:end-4])
   end

end

####################################################
filename = "box*outs";
dir = "/Users/hyzhou/Ganymede/run_mercury_80s_newbox/GM";

filename = "y*outs";
dir = "/Users/hyzhou/Ganymede/run_mercury_80s_newbox/GM";
