# Functions for the generation of a dipole field.
#
# Modified from the original Python version by
# Los Alamos National Security, LLC. 2010

using PyPlot

"""
	b_mag(x, y)

For a position `x`, `y` in units of planetary radius, return the strength of the
 dipole magnetic field in nT.
"""
function b_mag(x::AbstractFloat, y::AbstractFloat)
	r = @. √(x^2 + y^2)
	cosy = y./r
	B = @. 30400. * sqrt(1+3*cosy^2)/r^3
end

"""
	b_hat(x, y)
For given parameters, return two arrays, `x` and `y`, corresponding to the x and
 y components of b_hat for a dipole field. The grid is organized in meshgrid
(default) or ndgrid format.
"""
function b_hat(x, y; gridType="meshgrid")
   if gridType == "meshgrid"
	   xgrid = [i for j in y, i in x]
      ygrid = [j for j in y, i in x]
   else
      xgrid = [i for i in x, j in y]
      ygrid = [j for i in x, j in y]
   end

	r = @. sqrt(xgrid^2 + ygrid^2)
	cosy = ygrid./r
	sinx = xgrid./r

	denom = @. sqrt(1.0 + 3.0*cosy^2)

	br = @. 2.0 * cosy / denom
	bθ =          sinx ./ denom

	bx = @. br*sinx + bθ*cosy
	by = @. br*cosy - bθ*sinx

	return bx, by
end

"""
	b_line(x, y; npoint=30)

For a starting X, Y point return x and y vectors that trace the dipole field
line that passes through the given point.
"""
function b_line(x::AbstractFloat, y::AbstractFloat; npoints=30)
	r = sqrt(x^2 + y^2)
	try
		theta = atan(x/y)
	catch
		@warn "ZeroDivisionError"
		theta = pi/2.0
	end

	R = r/(sin(theta)^2)

	if x < 0
		theta = π:π/npoints:2.0π
	else
		theta = 0:π/npoints:π
	end

	r_vec = @. R * sin(theta)^2
	x_out = @. r_vec * sin(theta)
	y_out = @. r_vec * cos(theta)

	return x_out, y_out
end

"""A quick test of the dipole field functions."""
function test_dipole()

   x = -100.0:5.0:101.0
   y = -100.0:5.0:101.0

   x_vec, y_vec = b_hat(x,y)

   fig = plt.figure(figsize=(10,8))
   ax1 = plt.subplot(111)

   ax1.quiver(x, y, x_vec, y_vec)

   for i in -120:10:121
      x,y = b_line(float(i), 0.0, npoints=100)
      ax1.plot(x, y, "b")
   end
   for theta in π/2.0:π/100.0:3.0π/2.0
      x = sin(theta)
      y = cos(theta)
      x,y = b_line(x, y, npoints=100)
      ax1.plot(x, y, "r")
   end
   ax1.set_xlim([-100, 100])
   ax1.set_ylim([-100, 100])
   plt.title("Unit vectors for an arbitrary dipole field")

   return true
end
