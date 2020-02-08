# Field tracing routines for a streamline through a 2D vector field.
# Original version in C from LANL. Reimplement in Julia.
#
# Usage:
# 1. Run as it is in pure Julia.
# 2. Compile the standalone C code ctrace2d.c into a dynamic library and add to
#    path.
#
# Calling the native functions in Julia is about 5 times faster than calling the
# dynmic C library.
#
# Modified from [SpacePy](https://github.com/spacepy/spacepy)
# Hongyang Zhou, hyzhou@umich.edu 01/26/2020

using PyPlot, Random

export trace2d, trace2d_rk4, trace2d_eul, trace3d, trace3d_eul,
   test_trace_asymptote, test_trace_dipole, select_seeds

include("dipole.jl")

"""
	bilin_reg(x, y, Q00, Q01, Q10, Q11)

Bilinear interpolation for x1,y1=(0,0) and x2,y2=(1,1)
Q's are surrounding points such that Q00 = F[0,0], Q10 = F[1,0], etc.
"""
function bilin_reg(x, y, Q00, Q01, Q10, Q11)
   fout =
      Q00*(1.0-x)*(1.0-y) +
      Q10* x *    (1.0-y) +
      Q01* y *    (1.0-x) +
      Q11* x * y
end

"""
	DoBreak(iloc, jloc, iSize, jSize)

Check to see if we should break out of an integration.
"""
function DoBreak(iloc::Int, jloc::Int, iSize::Int, jSize::Int)
   ibreak = false
   if iloc ≥ iSize-1 || jloc ≥ jSize-1; ibreak = true end
   if iloc < 0 || jloc < 0; ibreak = true end
   return ibreak
end

function DoBreak(iloc::Int, jloc::Int, kloc::Int, iSize::Int, jSize::Int,
   kSize::Int)
   ibreak = false
   if iloc ≥ iSize-1 || jloc ≥ jSize-1 || kloc ≥ kSize-1; ibreak = true end
   if iloc < 0 || jloc < 0 || kloc < 0; ibreak = true end
   return ibreak
end

"""Create unit vectors of field."""
function make_unit!(iSize::Int, jSize::Int, ux, uy)
   for i = 1:iSize*jSize
      magnitude = sqrt(ux[i]^2 + uy[i]^2)
      ux[i] /= magnitude
      uy[i] /= magnitude
   end
end

function make_unit!(iSize::Int, jSize::Int, kSize::Int, ux, uy, uz)
   for i = 1:iSize*jSize*kSize
      magnitude = sqrt(ux[i]^2 + uy[i]^2 + uz[i]^2)
      ux[i] /= magnitude
      uy[i] /= magnitude
      uz[i] /= magnitude
   end
end

"""
	grid_interp!(x, y, field, xloc, yloc, xsize, ysize)

Interpolate a value at (x,y) in a field. `xloc` and `yloc` are indexes for x,y
locations (zero-based). `xsize` and `ysize` are the sizes of field in X and Y.
"""
grid_interp!(x, y, field, xloc, yloc, xsize, ysize) =
   bilin_reg(x-xloc, y-yloc,
   field[yloc*xsize+xloc+1],
   field[(yloc+1)*xsize+xloc+1],
   field[yloc*xsize+xloc+2],
   field[(yloc+1)*xsize+xloc+2])

"""
	Euler!(iSize,jSize, maxstep, ds, xstart,ystart, xGrid,yGrid, ux,uy, x,y)
Simple 2D tracing using Euler's method. Super fast but not super accurate.
# Arguments
- `iSize::Int,jSize::Int`: grid size.
- `maxstep::Int`: max steps.
- `ds::Float64`: step size.
- `xstart::Float64, ystart::Float64`: starting location.
- `xGrid::Array{Float64,2},yGrid::Array{Float64,2}`: actual coord system.
- `ux::Array{Float64,2},uy::Array{Float64,2}`: field to trace through.
- `x::Vector{Float64},y::Vector{Float64}`: x, y of result stream.
"""
function Euler!(iSize::Int, jSize::Int, maxstep::Int, ds,
   xstart, ystart, xGrid, yGrid, ux, uy, x, y)

   # Get starting points in normalized/array coordinates
   dx = xGrid[2] - xGrid[1]
   dy = yGrid[2] - yGrid[1]
   x[1] = (xstart-xGrid[1]) / dx
   y[1] = (ystart-yGrid[1]) / dy

   # Create unit vectors from full vector field
   make_unit!(iSize, jSize, ux, uy)

   nstep = 0
   # Perform tracing using Euler's method
   for n = 1:maxstep-1
      # Find surrounding points
      xloc = floor(Int, x[n])
      yloc = floor(Int, y[n])

      # Break if we leave the domain
      if DoBreak(xloc, yloc, iSize, jSize)
         nstep = n; break
      end

      # Interpolate unit vectors to current location
      fx = grid_interp!(x[n], y[n], ux, xloc, yloc, iSize, jSize)
      fy = grid_interp!(x[n], y[n], uy, xloc, yloc, iSize, jSize)

      # Detect NaNs in function values
      if isnan(fx) || isnan(fy) || isinf(fx) || isinf(fy)
         nstep = n
         break
      end

      # Perform single step
      x[n+1] = x[n] + ds * fx
      y[n+1] = y[n] + ds * fy

      nstep = n
   end

   # Return traced points to original coordinate system.
   for i = 1:nstep
      x[i] = x[i]*dx + xGrid[1]
      y[i] = y[i]*dy + yGrid[1]
   end
   return nstep
end


"""
	Euler!(iSize, jSize, kSize, maxstep, ds, xstart, ystart, zstart,
      xGrid, yGrid, zGrid ux, uy, uz, x, y, z)
Simple 3D tracing using Euler's method.
# Arguments
- `iSize::Int,jSize::Int,kSize::Int`: grid size.
- `maxstep::Int`: max steps.
- `ds::Float64`: step size.
- `xstart::Float64, ystart::Float64, zstart::Float64`: starting location.
- `xGrid::Array{Float64,2},yGrid::Array{Float64,2},yGrid::Array{Float64,2}`: actual coord system.
- `ux::Array{Float64,2},uy::Array{Float64,2},uz::Array{Float64,2}`: field to trace through.
- `x::Vector{Float64},y::Vector{Float64},z::Vector{Float64}`: x, y, z of result stream.
"""
function Euler!(iSize::Int, jSize::Int, kSize::Int, maxstep::Int, ds,
   xstart, ystart, zstart, xGrid, yGrid, zGrid, ux, uy, uz, x, y, z)

   # Get starting points in normalized/array coordinates
   dx = xGrid[2] - xGrid[1]
   dy = yGrid[2] - yGrid[1]
   dz = zGrid[2] - zGrid[1]
   x[1] = (xstart-xGrid[1]) / dx
   y[1] = (ystart-yGrid[1]) / dy
   z[1] = (zstart-zGrid[1]) / dz

   # Create unit vectors from full vector field
   make_unit!(iSize, jSize, kSize, ux, uy, uz)

   nstep = 0
   # Perform tracing using Euler's method
   for n = 1:maxstep-1
      # Find surrounding points
      xloc = floor(Int, x[n])
      yloc = floor(Int, y[n])
      zloc = floor(Int, z[n])

      # Break if we leave the domain
      if DoBreak(xloc, yloc, zloc, iSize, jSize, kSize)
         nstep = n; break
      end

      # Interpolate unit vectors to current location
      fx = grid_interp!(x[n], y[n], z[n], ux, xloc,yloc,zloc, iSize,jSize,kSize)
      fy = grid_interp!(x[n], y[n], z[n], uy, xloc,yloc,zloc, iSize,jSize,kSize)
      fy = grid_interp!(x[n], y[n], z[n], uz, xloc,yloc,zloc, iSize,jSize,kSize)

      # Detect NaNs in function values
      if any(isnan,[fx, fy, fz]) || any(isinf, [fx, fy, fz])
         nstep = n
         break
      end

      # Perform single step
      x[n+1] = x[n] + ds * fx
      y[n+1] = y[n] + ds * fy
      z[n+1] = z[n] + ds * fz

      nstep = n
   end

   # Return traced points to original coordinate system.
   for i = 1:nstep
      x[i] = x[i]*dx + xGrid[1]
      y[i] = y[i]*dy + yGrid[1]
      z[i] = z[i]*dz + zGrid[1]
   end
   return nstep
end


"""
	RK4!(iSize,jSize, maxstep, ds, xstart,ystart, xGrid,yGrid, ux,uy, x,y)

Fast and reasonably accurate 2D tracing with 4th order Runge-Kutta method and
constant step size `ds`.
"""
function RK4!(iSize::Int, jSize::Int, maxstep::Int, ds,
   xstart, ystart, xGrid, yGrid, ux, uy, x, y)

   # Get starting points in normalized/array coordinates
   dx = xGrid[2] - xGrid[1]
   dy = yGrid[2] - yGrid[1]
   x[1] = (xstart-xGrid[1]) / dx
   y[1] = (ystart-yGrid[1]) / dy

   # Create unit vectors from full vector field
   make_unit!(iSize, jSize, ux, uy)

   nstep = 0
   # Perform tracing using RK4
   for n = 1:maxstep-1
      # See Euler's method for more descriptive comments.
      # SUBSTEP #1
      xloc = floor(Int, x[n])
      yloc = floor(Int, y[n])
      if DoBreak(xloc, yloc, iSize, jSize); nstep = n; break end

      f1x = grid_interp!(x[n], y[n], ux, xloc, yloc, iSize, jSize)
      f1y = grid_interp!(x[n], y[n], uy, xloc, yloc, iSize, jSize)
      if isnan(f1x) || isnan(f1y) || isinf(f1x) || isinf(f1y)
         nstep = n; break
      end
      # SUBSTEP #2
      xpos = x[n] + f1x*ds/2.0
      ypos = y[n] + f1y*ds/2.0
      xloc = floor(Int, xpos)
      yloc = floor(Int, ypos)
      if DoBreak(xloc, yloc, iSize, jSize); nstep = n; break end

      f2x = grid_interp!(xpos, ypos, ux, xloc, yloc, iSize, jSize)
      f2y = grid_interp!(xpos, ypos, uy, xloc, yloc, iSize, jSize)

      if isnan(f2x) || isnan(f2y) || isinf(f2x) || isinf(f2y)
         nstep = n; break
      end
      # SUBSTEP #3
      xpos = x[n] + f2x*ds/2.0
      ypos = y[n] + f2y*ds/2.0
      xloc = floor(Int, xpos)
      yloc = floor(Int, ypos)
      if DoBreak(xloc, yloc, iSize, jSize); nstep = n; break end

      f3x = grid_interp!(xpos, ypos, ux, xloc, yloc, iSize, jSize)
      f3y = grid_interp!(xpos, ypos, uy, xloc, yloc, iSize, jSize)
      if isnan(f3x) || isnan(f3y) || isinf(f3x) || isinf(f3y)
         nstep = n; break
      end

      # SUBSTEP #4
      xpos = x[n] + f3x*ds
      ypos = y[n] + f3y*ds
      xloc = floor(Int, xpos)
      yloc = floor(Int, ypos)
      if DoBreak(xloc, yloc, iSize, jSize); nstep = n; break end

      f4x = grid_interp!(xpos, ypos, ux, xloc, yloc, iSize, jSize)
      f4y = grid_interp!(xpos, ypos, uy, xloc, yloc, iSize, jSize)
      if isnan(f4x) || isnan(f4y) || isinf(f4x) || isinf(f4y)
         nstep = n; break
      end

      # Peform the full step using all substeps
      x[n+1] = x[n] + ds/6.0 * (f1x + f2x*2.0 + f3x*2.0 + f4x)
      y[n+1] = y[n] + ds/6.0 * (f1y + f2y*2.0 + f3y*2.0 + f4y)

      nstep = n
   end

   # Return traced points to original coordinate system.
   for i = 1:nstep
      x[i] = x[i]*dx + xGrid[1]
      y[i] = y[i]*dy + yGrid[1]
   end
   return nstep
end


"""
	RK4!(iSize,jSize,kSize, maxstep, ds, xstart,ystart,zstart, xGrid,yGrid,zGrid,
      ux,uy,uz, x,y,z)

Fast and reasonably accurate 3D tracing with 4th order Runge-Kutta method and
constant step size `ds`.
"""
function RK4!(iSize::Int, jSize::Int, kSize::Int, maxstep::Int, ds,
   xstart, ystart, zstart, xGrid, yGrid, zGrid, ux, uy, uz, x, y, z)

   # Get starting points in normalized/array coordinates
   dx = xGrid[2] - xGrid[1]
   dy = yGrid[2] - yGrid[1]
   dz = zGrid[2] - zGrid[1]
   x[1] = (xstart-xGrid[1]) / dx
   y[1] = (ystart-yGrid[1]) / dy
   z[1] = (zstart-zGrid[1]) / dz

   # Create unit vectors from full vector field
   make_unit!(iSize, jSize, kSize, ux, uy, uz)

   nstep = 0
   # Perform tracing using RK4
   for n = 1:maxstep-1
      # See Euler's method for more descriptive comments.
      # SUBSTEP #1
      xloc = floor(Int, x[n])
      yloc = floor(Int, y[n])
      zloc = floor(Int, z[n])
      if DoBreak(xloc, yloc, zloc, iSize, jSize, kSize); nstep = n; break end

      f1x = grid_interp!(x[n],y[n],z[n], ux, xloc,yloc,zloc, iSize,jSize,kSize)
      f1y = grid_interp!(x[n],y[n],z[n], uy, xloc,yloc,zloc, iSize,jSize,kSize)
      f1z = grid_interp!(x[n],y[n],z[n], uz, xloc,yloc,zloc, iSize,jSize,kSize)
      if any(isnan,[f1x, f1y, f1z]) || any(isinf, [f1x, f1y, f1z])
         nstep = n; break
      end
      # SUBSTEP #2
      xpos = x[n] + f1x*ds/2.0
      ypos = y[n] + f1y*ds/2.0
      zpos = z[n] + f1z*ds/2.0
      xloc = floor(Int, xpos)
      yloc = floor(Int, ypos)
      zloc = floor(Int, zpos)
      if DoBreak(xloc, yloc, zloc, iSize, jSize, kSize); nstep = n; break end

      f2x = grid_interp!(xpos,ypos,zpos, ux, xloc,yloc,zloc, iSize,jSize,kSize)
      f2y = grid_interp!(xpos,ypos,zpos, uy, xloc,yloc,zloc, iSize,jSize,kSize)
      f2z = grid_interp!(xpos,ypos,zpos, uz, xloc,yloc,zloc, iSize,jSize,kSize)
      if any(isnan,[f2x, f2y, f2z]) || any(isinf, [f2x, f2y, f2z])
         nstep = n; break
      end
      # SUBSTEP #3
      xpos = x[n] + f2x*ds/2.0
      ypos = y[n] + f2y*ds/2.0
      zpos = z[n] + f2z*ds/2.0
      xloc = floor(Int, xpos)
      yloc = floor(Int, ypos)
      zloc = floor(Int, zpos)
      if DoBreak(xloc, yloc, zloc, iSize, jSize, kSize); nstep = n; break end

      f3x = grid_interp!(xpos,ypos,zpos, ux, xloc,yloc,zloc, iSize,jSize,kSize)
      f3y = grid_interp!(xpos,ypos,zpos, uy, xloc,yloc,zloc, iSize,jSize,kSize)
      f3z = grid_interp!(xpos,ypos,zpos, uz, xloc,yloc,zloc, iSize,jSize,kSize)
      if any(isnan,[f3x, f3y, f3z]) || any(isinf, [f3x, f3y, f3z])
         nstep = n; break
      end

      # SUBSTEP #4
      xpos = x[n] + f3x*ds
      ypos = y[n] + f3y*ds
      zpos = z[n] + f3z*ds
      xloc = floor(Int, xpos)
      yloc = floor(Int, ypos)
      zloc = floor(Int, zpos)
      if DoBreak(xloc, yloc, zloc, iSize, jSize, kSize); nstep = n; break end

      f4x = grid_interp!(xpos,ypos,zpos, ux, xloc,yloc,zloc, iSize,jSize,kSize)
      f4y = grid_interp!(xpos,ypos,zpos, uy, xloc,yloc,zloc, iSize,jSize,kSize)
      f4z = grid_interp!(xpos,ypos,zpos, uz, xloc,yloc,zloc, iSize,jSize,kSize)
      if any(isnan,[f4x, f4y, f4z]) || any(isinf, [f4x, f4y, f4z])
         nstep = n; break
      end

      # Peform the full step using all substeps
      x[n+1] = x[n] + ds/6.0 * (f1x + f2x*2.0 + f3x*2.0 + f4x)
      y[n+1] = y[n] + ds/6.0 * (f1y + f2y*2.0 + f3y*2.0 + f4y)
      z[n+1] = z[n] + ds/6.0 * (f1z + f2z*2.0 + f3z*2.0 + f4z)

      nstep = n
   end

   # Return traced points to original coordinate system.
   for i = 1:nstep
      x[i] = x[i]*dx + xGrid[1]
      y[i] = y[i]*dy + yGrid[1]
      z[i] = z[i]*dz + zGrid[1]
   end
   return nstep
end


"""
	trace2d_rk4(fieldx, fieldy, xstart, ystart, gridx, gridy;
		maxstep=20000, ds=0.01, gridType="meshgrid", direction="both")

Given a 2D vector field, trace a streamline from a given point to the edge of
the vector field. The field is integrated using Runge Kutta 4. Slower than
Euler, but more accurate. The higher accuracy allows for larger step sizes `ds`.
 Only valid for regular grid with coordinates `gridx`, `gridy`. If `gridx` and
`gridy` are not given, assume that `xstart` and `ystart` are normalized
coordinates (e.g., position in terms of array indices.???)
The field can be in both `meshgrid` (default) or `ndgrid` format.
Supporting `direction` for {"both","forward","backward"}.
"""
function trace2d_rk4(fieldx, fieldy, xstart, ystart, gridx, gridy;
   maxstep=20000, ds=0.01, gridType="meshgrid", direction="both")

   xt = Vector{eltype(fieldx)}(undef,maxstep) # output x
   yt = Vector{eltype(fieldy)}(undef,maxstep) # output y

   nx, ny = size(gridx)[1], size(gridy)[1]

   gx, gy = collect(gridx), collect(gridy)

   if gridType == "ndgrid" # ndgrid
      fx = fieldx
      fy = fieldy
   else # meshgrid
      fx = permutedims(fieldx)
      fy = permutedims(fieldy)
   end

   #=
   if eltype(fieldx) == Float64
      npoints = ccall((:cRk4,"libtrace.so"),Cint,
      (Cint,Cint,Cint,Float64,Float64,Float64,Ptr{Cdouble},Ptr{Cdouble},
      Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
      nx, ny, maxstep, ds, xstart, ystart, gx, gy, fx, fy, xt, yt)
   elseif eltype(fieldx) == Float32
      npoints = ccall((:cRk4,"libtrace_single.so"),Cint,
      (Cint,Cint,Cint,Float32,Float32,Float32,Ptr{Cfloat},Ptr{Cfloat},
      Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}),
      Int32(nx), Int32(ny), Int32(maxstep), Float32(ds), xstart, ystart, gx, gy,
      fx, fy, xt, yt)
   end
   =#

   if direction == "forward"
      npoints = RK4!(nx,ny, maxstep, ds, xstart,ystart, gx,gy, fx,fy, xt,yt)
   elseif direction == "backward"
      npoints = RK4!(nx,ny, maxstep, ds, xstart,ystart, gx,gy, -fx,-fy, xt,yt)
   else
      n1 = RK4!(nx,ny, maxstep, ds, xstart,ystart, gx,gy,-fx,-fy, xt,yt)
      xt[n1:-1:1] = xt[1:n1]
      yt[n1:-1:1] = yt[1:n1]

      x2 = Vector{eltype(fieldx)}(undef,maxstep-n1)
      y2 = Vector{eltype(fieldx)}(undef,maxstep-n1)
      n2 = RK4!(nx,ny, maxstep-n1, ds, xstart,ystart, gx,gy, fx,fy, x2,y2)
      xt[n1+1:n1+n2-1] = x2[2:n2]
      yt[n1+1:n1+n2-1] = y2[2:n2]
      npoints = n1 + n2 - 1
   end

   return xt[1:npoints], yt[1:npoints]
end

"""
	trace2d_eul(fieldx, fieldy, xstart, ystart, gridx, gridy;
		maxstep=20000, ds=0.01, gridType="meshgrid", direction="both")

Given a 2D vector field, trace a streamline from a given point to the edge of
the vector field. The field is integrated using Euler's method. While this is
faster than rk4, it is less accurate. Only valid for regular grid with
coordinates `gridx`, `gridy`. If gridx and gridy are not given, assume that
`xstart` and `ystart` are normalized coordinates (e.g. position in terms of
array indices.)??? The field can be in both `meshgrid` (default) or `ndgrid`
format.
Supporting `direction` for {"both","forward","backward"}.
"""
function trace2d_eul(fieldx, fieldy, xstart, ystart, gridx, gridy;
   maxstep=20000, ds=0.01, gridType="meshgrid", direction="both")

   xt = Vector{eltype(fieldx)}(undef,maxstep) # output x
   yt = Vector{eltype(fieldy)}(undef,maxstep) # output y

   nx, ny = size(gridx)[1], size(gridy)[1]

   gx, gy = collect(gridx), collect(gridy)

   if gridType == "ndgrid" # ndgrid
      fx = fieldx
      fy = fieldy
   else # meshgrid
      fx = permutedims(fieldx)
      fy = permutedims(fieldy)
   end

   #=
   if eltype(fieldx) == Float64
      npoints = ccall((:cEuler,"libtrace.so"),Cint,
      (Cint,Cint,Cint,Float64,Float64,Float64,Ptr{Cdouble},Ptr{Cdouble},
      Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
      nx, ny, maxstep, ds, xstart, ystart, gx, gy, fx, fy, xt, yt)
   elseif eltype(fieldx) == Float32
      npoints = ccall((:cEuler,"libtrace_single.so"),Cint,
      (Cint,Cint,Cint,Float32,Float32,Float32,Ptr{Cfloat},Ptr{Cfloat},
      Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}),
      Int32(nx), Int32(ny), Int32(maxstep), Float32(ds), xstart, ystart, gx, gy,
      fx, fy, xt, yt)
   end
   =#

   if direction == "forward"
      npoints = Euler!(nx,ny, maxstep, ds, xstart,ystart, gx,gy, fx,fy, xt,yt)
   elseif direction == "backward"
      npoints = Euler!(nx,ny, maxstep, ds, xstart,ystart, gx,gy, -fx,-fy, xt,yt)
   else
      n1 = Euler!(nx,ny, maxstep, ds, xstart,ystart, gx,gy,-fx,-fy, xt,yt)
      xt[n1:-1:1] = xt[1:n1]
      yt[n1:-1:1] = yt[1:n1]

      x2 = Vector{eltype(fieldx)}(undef,maxstep-n1)
      y2 = Vector{eltype(fieldx)}(undef,maxstep-n1)
      n2 = Euler!(nx,ny, maxstep-n1, ds, xstart,ystart, gx,gy, fx,fy, x2,y2)
      xt[n1+1:n1+n2-1] = x2[2:n2]
      yt[n1+1:n1+n2-1] = y2[2:n2]
      npoints = n1 + n2 - 1
   end

   return xt[1:npoints], yt[1:npoints]
end

trace2d(fieldx, fieldy, xstart, ystart, gridx, gridy) =
   trace2d_rk4(fieldx, fieldy, xstart, ystart, gridx, gridy)

"""
	trace3d_eul(fieldx, fieldy, fieldz, xstart, ystart, zstart, gridx, gridy,
      gridz; maxstep=20000, ds=0.01)

Given a 3D vector field, trace a streamline from a given point to the edge of
the vector field. The field is integrated using Euler's method. Only valid for
regular grid with coordinates `gridx`, `gridy`, `gridz`.
The field can be in both `meshgrid` (default) or `ndgrid` format.
"""
function trace3d_eul(fieldx, fieldy, fieldz, xstart, ystart, zstart, gridx,
   gridy, gridz; maxstep=20000, ds=0.01, gridType="meshgrid")

   xt = Vector{eltype(fieldx)}(undef,maxstep) # output x
   yt = Vector{eltype(fieldy)}(undef,maxstep) # output y
   zt = Vector{eltype(fieldz)}(undef,maxstep) # output z

   nx, ny, nz = size(gridx)[1], size(gridy)[1], size(gridz)[1]

   gx, gy, gz = collect(gridx), collect(gridy), collect(gridz)

   if gridType == "ndgrid" # ndgrid
      fx, fy, fz = fieldx, fieldy, fieldz
   else # meshgrid
      fx, fy, fz = permutedims(fieldx), permutedims(fieldy), permutedims(fieldz)
   end

   npoints = Euler!(nx,ny,nz, maxstep, ds, xstart,ystart,zstart, gx,gy,gz,
      fx,fy,fz, xt,yt,zt)

   return xt[1:npoints], yt[1:npoints], zt[1:npoints]
end

trace3d(fieldx, fieldy, fieldz, xstart, ystart, zstart, gridx, gridy, gridz) =
   trace3d_eul(fieldx, fieldy, fieldz, xstart, ystart, zstart,
   gridx, gridy, gridz)

"""
	select_seeds(x, y, nSeed=100)

Generate `nSeed` seeding points randomly in the grid range. If you specify
`nSeed`, use the keyword input, otherwise it will be overloaded by the 3D
version seed generation.
"""
function select_seeds(x, y; nSeed=100)
   xmin,xmax = extrema(x)
   ymin,ymax = extrema(y)

   xstart = rand(MersenneTwister(0),nSeed)*(xmax-xmin) .+ xmin
   ystart = rand(MersenneTwister(1),nSeed)*(ymax-ymin) .+ ymin
   seeds = zeros(eltype(x[1]),2,nSeed)
   for i = 1:length(xstart)
      seeds[1,i] = xstart[i]
      seeds[2,i] = ystart[i]
   end
   return seeds
end

function select_seeds(x, y, z; nSeed=100)
   xmin,xmax = extrema(x)
   ymin,ymax = extrema(y)
   zmin,zmax = extrema(y)

   xstart = rand(MersenneTwister(0),nSeed)*(xmax-xmin) .+ xmin
   ystart = rand(MersenneTwister(1),nSeed)*(ymax-ymin) .+ ymin
   zstart = rand(MersenneTwister(2),nSeed)*(zmax-zmin) .+ zmin
   seeds = zeros(eltype(x[1]),3,length(xstart))
   for i = 1:length(xstart)
      seeds[1,i] = xstart[i]
      seeds[2,i] = ystart[i]
      seeds[3,i] = zstart[i]
   end
   return seeds
end

"""
   test_trace_asymptote(IsSingle)

Test streamline tracing by plotting vectors and associated streamlines
through a simple velocity field where Vx=x, Vy=-y.
Support test for single and double precision.
"""
function test_trace_asymptote(IsSingle=false)

   # Start by creating a velocity vector field.
   if IsSingle
      xmax, ymax = 200.0f0, 20.0f0
      x = -10.0f0:0.25f0:xmax
      y = -10.0f0:0.25f0:ymax
   else
      xmax, ymax = 200.0, 20.0
      x = -10.0:0.25:xmax
      y = -10.0:0.25:ymax
   end

   xgrid = [i for j in y, i in x]
   ygrid = [j for j in y, i in x]

   if IsSingle
      vx = xgrid * 1.0f0
      vy = ygrid * -1.0f0

      xstart = 1.0f0
      ystart = 10.0f0
   else
      vx = xgrid * 1.0
      vy = ygrid * -1.0

      xstart = 1.0
      ystart = 10.0
   end

   x1, y1 = trace2d_rk4(vx, vy, xstart, ystart, x, y, ds=0.1)
   x2, y2 = trace2d_rk4(vx, vy, xstart, ystart, x, y, ds=0.5)
   x3, y3 = trace2d_rk4(vx, vy, xstart, ystart, x, y, ds=1.0)
   x4, y4 = trace2d_eul(vx, vy, xstart, ystart, x, y, ds=0.1)
   x5, y5 = trace2d_eul(vx, vy, xstart, ystart, x, y, ds=0.5)
   x6, y6 = trace2d_eul(vx, vy, xstart, ystart, x, y, ds=1.0)

   # analytical solution const = x*y
   c = xstart * ystart
   if IsSingle
      x_anly = 1.0f0:0.001f0:xmax
   else
      x_anly = 1.0:0.001:xmax
   end
   y_anly = c ./ x_anly

   fig = plt.figure(figsize=(8,6))
   ax1 = plt.subplot(111)
   ax1.plot(x_anly, y_anly, "r--", label="Analytic",linewidth=3.0)
   ax1.plot(x1, y1, "b",   label="RK4 ds=0.1", linewidth=1.5)
   ax1.plot(x2, y2, "b--", label="RK4 ds=0.5", linewidth=1.5)
   ax1.plot(x3, y3, "b:",  label="RK4 ds=1.0", linewidth=1.5)
   ax1.plot(x4, y4, "g", label="Euler ds=0.1", linewidth=.75)
   ax1.plot(x5, y5, "g--", label="Euler ds=0.5", linewidth=.75)
   ax1.plot(x6, y6, "g:",  label="Euler ds=1.0", linewidth=.75)
   ax1.legend(loc="upper right")
   if IsSingle
      ax1.set_title("Runge Kutta 4 vs Euler's: Asymptotic Field, Single Precision")
   else
      ax1.set_title("Runge Kutta 4 vs Euler's: Asymptotic Field, Double Precision")
   end
   ax1.set_xlabel("Normalized X Coordinate")
   ax1.set_ylabel("Normalized Y Coordinate")
   ax1.set_xlim([3, 30])
   ax1.set_ylim([0.25, 2.5 ])

   # Annotate plot.
   ax1.annotate("Euler's method diverges strongly\nalong curves except when"*
      " taking \nvery small steps.  RK4 is\nfar more accurate "*
      "for all dS.",
      xy=(10.5,0.83), xycoords="data", xytext=(9,.3),
      arrowprops=Dict(:fc=>"black",:shrink=>0.05),
      horizontalalignment="center")
   ax1.annotate("Tracing begins at x=1, y=10.",
      xy=(4.8,2.45), xycoords="data", xytext=(6,1.9),
      arrowprops=Dict(:fc=>"black",:shrink=>0.05))
   ax1.annotate("The trace continues until \n"*
      "x=200.  At that point, Euler \n"*
      "dS=0.1 and RK4 dS=1.0 converge\n"*
      " at the same point, despite \n"*
      "the 10X difference in step size.",
      xy=(29,.5), xycoords="data", xytext=(20,1),
      arrowprops=Dict(:fc=>"black",:shrink=>0.005),
      horizontalalignment="center")

   return true
end

"""Trace field lines through a dipole field."""
function test_trace_dipole()
   # Now do dipole magnetic field.
   # Start by creating a field of unit vectors...
   x = -100.0:5.0:101.0
   y = -100.0:5.0:101.0

   bx, by = b_hat(x, y)

   # New figure
   fig2 = plt.figure(figsize=(12,6))
   fig2.subplots_adjust(wspace=0.15, left=0.08, right=0.94)
   ax2 = plt.subplot(121)
   ax3 = plt.subplot(322)
   ax4a= plt.subplot(347)
   ax4b= plt.subplot(348)
   ax5 = plt.subplot(326)
   ax2.quiver(x,y, bx, by, units="x", pivot="middle")
   ax3.quiver(x,y, bx, by, units="x", pivot="tip")
   ax4a.quiver(x,y,bx, by, units="x", pivot="tip")
   ax4b.quiver(x,y,bx, by, units="x", pivot="tip")
   ax5.quiver(x,y, bx, by, units="x", pivot="tip")

   # Trace through this field.
   xstart = 10.0 # ystart = 25.0
   ds = 0.1
   for ystart in 0.:5.:31.
      x1, y1 = trace2d_rk4(bx, by, xstart, ystart, x, y, ds=ds)
      l1 = ax2.plot(x1,y1,"b")[1]
      ax3.plot(x1,y1,"b"); ax4b.plot(x1,y1,"b")
      ax5.plot(x1,y1,"b"); ax4a.plot(x1,y1,"b")
      x2, y2 = trace2d_eul(bx, by, xstart, ystart, x, y, ds=ds)
      l2 = ax2.plot(x2,y2,"r")[1]
      ax3.plot(x2,y2,"r"); ax4b.plot(x2,y2,"r")
      ax5.plot(x2,y2,"r"); ax4a.plot(x2,y2,"r")
      x3, y3 = b_line(xstart, ystart, npoints=300)
      l3 = ax2.plot(x3,y3,"k--")[1]
      ax3.plot(x3,y3,"k--"); ax4b.plot(x3,y3,"k--")
      ax5.plot(x3,y3,"k--"); ax4a.plot(x3,y3,"k--")

      ax2.set_xlim([-2,  100])
      ax2.set_ylim([-30, 100])
      ax2.set_title("Full View")
      ax2.set_xlabel("Normalized X Coordinate")
      ax2.set_ylabel("Normalized Y Coordinate")
      ax2.legend((l1, l2, l3),("RK4", "Euler", "Analytical"), loc="upper left")

      ax3.set_title("Zoomed Views")
      ax3.set_xlim([8.5, 17.5])
      ax3.set_ylim([3, 33])
      goodpos = ax3.get_position()

      ax4a.set_xlim([20,30])
      ax4a.set_ylim([-12,12])
      pos = ax4a.get_position()
      pos.x0 = goodpos.x0
      pos.x1 = pos.x0 + (goodpos.x1-goodpos.x0)/2.0 -0.01
      ax4a.set_position(pos)

      ax4b.set_xlim([50,60])
      ax4b.set_ylim([-12,12])
      pos = ax4b.get_position()
      pos.x0 = goodpos.x0 + (goodpos.x1-goodpos.x0)/2.0 +0.01
      ax4b.set_position(pos)
      ax4b.yaxis.set_ticklabels([])

      ax5.set_xlim([1,7])
      ax5.set_ylim([-7,-3])
      ax5.set_xlabel("Normalized X Coordinate")

      fig2.suptitle("RK4 vs Euler's Method: Dipole Field for dS="*
         string(round(ds,digits=3)))
   end
   return true
end
