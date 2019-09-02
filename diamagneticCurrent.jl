using PyCall, PyPlot#, ScatteredInterpolation

include("VisAna.jl")
using .VisAna
filename = "3d.dat"
filehead, data, filelist = readdata(filename,verbose=false);

function get_diamagnetic_current(filehead::Dict, data::Array{T,2}) where
   T<:AbstractFloat

   @views X = data[1,:]
   @views Y = data[2,:]
   @views Z = data[3,:]
   @views Bx0 = data[11,:]
   @views By0 = data[12,:]
   @views Bz0 = data[13,:]
   @views P0 = data[14,:]

   points = hcat(X,Y,Z);
   #points = [X Y Z]' # for ScatteredInterpolation

   np = pyimport("numpy")
   scipy = pyimport("scipy")
   interpolate = pyimport("scipy.interpolate")
   griddata = interpolate.griddata

   #=
   xMin = minimum(X); xMax = maximum(X)
   yMin = minimum(Y); yMax = maximum(Y)
   zMin = minimum(Z); zMax = maximum(Z)
   =#

   xMin, xMax = -3.0, 3.0
   yMin, yMax = -3.0, 3.0
   zMin, zMax = -3.0, 3.0

   nX, nY, nZ = 100, 3, 100

   s = pybuiltin(:slice)
   XYZ = get(np.mgrid,
      (s(xMin,xMax,nX*1im), s(yMin,yMax,nY*1im), s(zMin,zMax,nZ*1im)))
   grid_x, grid_y, grid_z = XYZ[1,:,:,:], XYZ[2,:,:,:], XYZ[3,:,:,:]


   #= # The major issue with ScatteredInterpolation is memory!
   XYZ = reshape(XYZ,3,nX*nY*nZ)

   # Let me try ScatteredInterpolation in Julia
   itp = interpolate(Shepard(), points, Bx0)
   Bx = evaluate(itp, XYZ)
   Bx = reshape(Bx,nX,nY,nZ)

   itp = interpolate(Shepard(), points, By0)
   By = evaluate(itp, XYZ)
   By = reshape(By,nX,nY,nZ)

   itp = interpolate(Shepard(), points, Bz0)
   Bz = evaluate(itp, XYZ)
   Bz = reshape(Bz,nX,nY,nZ)

   itp = interpolate(Shepard(), points, P0)
   P = evaluate(itp, XYZ)
   P = reshape(P,nX,nY,nZ)
   =#

   Bx = griddata(points, Bx0, (grid_x, grid_y, grid_z), method="linear")
   By = griddata(points, By0, (grid_x, grid_y, grid_z), method="linear")
   Bz = griddata(points, Bz0, (grid_x, grid_y, grid_z), method="linear")
   P  = griddata(points, P0,  (grid_x, grid_y, grid_z), method="linear")

   px, py, pz = np.gradient(P) # ▽P

   jDiaMag  = Array{Float64,4}(undef,3,nX,nY,nZ)
   jCurrent = Array{Float64,3}(undef,nX,nY,nZ)

   # j = B×▽p/B^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag = sqrt(Bx[i,j,k]^2 + Bx[i,j,k]^2 + Bz[i,j,k]^2)
      jDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) / Bmag
      jDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) / Bmag
      jDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) / Bmag
      jCurrent[i,j,k] =
      sqrt(jDiaMag[1,i,j,k]^2 + jDiaMag[2,i,j,k]^2 + jDiaMag[3,i,j,k]^2)
   end

   x = range(xMin, stop=xMax, length=nX)
   z = range(zMin, stop=zMax, length=nZ)


   contourf(x,z,jCurrent[:,2,:]',50)
   colorbar()
   xlabel("x")
   ylabel("z")

end

get_diamagnetic_current(filehead[1], data)
