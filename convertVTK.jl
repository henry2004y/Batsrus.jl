# Convert SWMF GM outputs to vtk files.
#
# Hongyang Zhou, hyzhou@umich.edu 07/23/2019

using WriteVTK
try
   using MATLAB
catch
   println("MATLAB not found path...")
end
include("VisAna.jl")
using .VisAna

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

mypath = "."
mykey  = "box_test.out"

filenames = searchdir(mypath,mykey)

function convertBox2VTK_matlab(filenames::Array{String,1})

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

      outname = filename[1:end-4]

      # Rectilinear or structured grid
      outfiles = vtk_grid(outname, x, y, z) do vtk
         vtk_point_data(vtk, P, "P")
         vtk_point_data(vtk, B, "B")
      end
   end

end

function convertBox2VTK(filenames::Array{String,1}, gridType::Int64=1)

   for filename in filenames
      filehead, data, filelist = readdata(filename,verbose=false)

      func = "bx"
      func_ = findfirst(x->x==func, lowercase.(filehead[1][:wnames]))
      if isnothing(func_) @error "Couldn't find variable $(func)!" end
      B = @view data[1].w[:,:,:,func_:func_+2]
      B = permutedims(B, [4,1,2,3])

      #=
      func = "P"
      func_ = findfirst(x->x==func, filehead[1][:wnames])
      if isnothing(func_) @error "Couldn't find variable $(func)!" end
      P = data[1].w[:,:,:,func_]
      =#

      outname = filename[1:end-4]

      if gridType == 1 # rectilinear grid
         x = @view data[1].x[:,1,1,1]
         y = @view data[1].x[1,:,1,2]
         z = @view data[1].x[1,1,:,3]

         outfiles = vtk_grid(outname, x,y,z) do vtk
            #vtk_point_data(vtk, P, "P")
            vtk_point_data(vtk, B, "B")
         end
      elseif gridType == 2 # structured grid
         xyz = permutedims(data[1].x, [4,1,2,3])

         outfiles = vtk_grid(outname, xyz) do vtk
            #vtk_point_data(vtk, P, "P")
            vtk_point_data(vtk, B, "B")
         end
      elseif gridType == 3 # unstructured grid, not finished
         vtkfile = vtk_grid(outname, points, cells)
         vtk_cell_data(vtkfile, cdata, "my_cell_data")
         outfiles = vtk_save(vtkfile)
      end
      println(filename," finished conversion.")
   end

end

function convertTec2VTK()
   filename = "3d.dat"
   head, data, connectivity  = readtecdata(filename,false)

   points = @view data[1:3,:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:nDim] == 3
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:nDim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   # 777MB in Ascii to 147MB in VTK binary;
   # What about preplot? 234MB
   vtkfile = vtk_grid("test_unstructured", points, cells)

   rho = @view data[4,:]
   p = @view data[14,:]

   vtk_point_data(vtkfile, rho, "Rho")
   vtk_point_data(vtkfile, p, "P")

   outfiles = vtk_save(vtkfile)
end

function test_bin()
   filename = "3d_mhd.dat"

   head, data, connectivity  = readtecdata(filename,true)

   points = @view data[1:3,:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:nDim] == 3
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:nDim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid("test_unstructured", points, cells)

   rho = @view data[4,:]
   p = @view data[14,:]

   vtk_point_data(vtkfile, rho, "Rho")
   vtk_point_data(vtkfile, p, "P")

   outfiles = vtk_save(vtkfile)

end