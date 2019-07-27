# Convert SWMF GM outputs to vtk files.
# Right now only works for Cartesian grids.
#
# Hongyang Zhou, hyzhou@umich.edu 07/23/2019

using WriteVTK, MATLAB

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

mypath = "."
mykey  = "box"

filenames = searchdir(mypath,mykey)

function convertBox2VTK(filenames::Array{String,1})

   for filename in filenames
      filehead, data = mxcall(:read_data,2,filename)

      data = data["file1"]

      x = data["x"][:,:,:,1]
      y = data["x"][:,:,:,2]
      z = data["x"][:,:,:,3]

      func = "Bx"
      func_ = findfirst(x->x==func, filehead["wnames"])
      if isnothing(func_) @error "Couldn't find variable $(func)!" end
      B = data["w"][:,:,:,func_:func_+2]
      B = permutedims(B, [4,1,2,3])

      func = "P"
      func_ = findfirst(x->x==func, filehead["wnames"])
      if isnothing(func_) @error "Couldn't find variable $(func)!" end
      P = data["w"][:,:,:,func_]

      outname = filename[1:19]

      # Rectilinear or structured grid
      outfiles = vtk_grid(outname, x, y, z) do vtk
         vtk_point_data(vtk, P, "P")
         vtk_point_data(vtk, B, "B")
      end
   end

end
