# Test of using user recipes from Plots.
#
# Hongyang Zhou, hyzhou@umich.edu 12/01/2019

using Plots

#@userplot mycontour

# Build a recipe which acts on a custom type.
# Notice that the function apply_recipe is returned.
# The recipe macro is just a convenience to build apply_recipe definitions.
@recipe function plotdata(data::Data, filehead::Dict, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)

	ndim = filehead[:ndim]

   x,w = data.x, data.w
   VarIndex_ = findfirst(x->x==var,filehead[:wnames])

   if ndim == 1
      y = w[:,VarIndex_]
      @series begin
         seriestype --> :path
         x, y
      end
   elseif ndim == 2
      if filehead[:gencoord] # Generalized coordinates
         X = vec(x[:,:,1])
         Y = vec(x[:,:,2])
         W = vec(w[:,:,VarIndex_])

         if any(abs.(plotrange) .== Inf)
            if plotrange[1] == -Inf plotrange[1] = minimum(X) end
            if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
            if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
            if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
         end

         # Create grid values first.
         xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
         yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
         # Perform linear interpolation of the data (x,y) on grid(xi,yi)
         # This part is not working because I don't know how to import
         # matplotlib!
         triang = matplotlib.tri.Triangulation(X,Y)
         interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
         Xi = [y for x in xi, y in yi]
         Yi = [x for x in xi, y in yi]
         wi = interpolator(Xi, Yi)
      else # Cartesian coordinates
         if all(isinf.(plotrange))
            xi = x[:,1,1]
            yi = x[1,:,2]
            wi = w[:,:,VarIndex_]
         else
            X = x[:,1,1]
            Y = x[1,:,2]
            if plotrange[1] == -Inf plotrange[1] = minimum(X) end
            if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
            if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
            if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

            W = w[:,:,VarIndex_]

            xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
            yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

            spline = Spline2D(X, Y, W)
            Xi = [i for i in xi, j in yi]
            Yi = [j for i in xi, j in yi]
            wi = spline(Xi[:], Yi[:])
            wi = reshape(wi, size(Xi))'
         end
      end

      # The GR backend contour plot only accept vector x,y!

      @series begin
         seriestype := :contour  # use --> if you don't want to force it
         xi, yi, wi
      end
   end
end

# Type annotations on keyword arguments not currently supported in recipes.
# Type information has been discarded
#=
@recipe function contourf(data::Data, filehead::Dict, var::String,
	levels::Int=0,
   plotrange::Vector{Float64}=[-Inf,Inf,-Inf,Inf], plotinterval::Float64=0.1)

   xi, yi, wi = getdata(data, filehead, var, plotrange, plotinterval)

	return xi, yi, wi
end
=#

"""Prepare data for passing to plotting functions."""
function getdata(data, filehead, var, plotrange, plotinterval)
   x,w = data.x, data.w
   ndim = filehead[:ndim]

   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(filehead[:wnames]))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

   if filehead[:gencoord] # Generalized coordinates
      X = vec(x[:,:,1])
      Y = vec(x[:,:,2])
      W = vec(w[:,:,VarIndex_])

      if any(abs.(plotrange) .== Inf)
         if plotrange[1] == -Inf plotrange[1] = minimum(X) end
         if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
         if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
         if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
      end

      # Create grid values first.
      xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
      yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
      # Perform linear interpolation of the data (x,y) on grid(xi,yi)
      triang = matplotlib.tri.Triangulation(X,Y)
      interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
      Xi = [y for x in xi, y in yi]
	   Yi = [x for x in xi, y in yi]
      wi = interpolator(Xi, Yi)
   else # Cartesian coordinates
      if all(isinf.(plotrange))
         xi = x[:,:,1]
         yi = x[:,:,2]
         wi = w[:,:,VarIndex_]
      else
		   X = x[:,1,1]
		   Y = x[1,:,2]
         if plotrange[1] == -Inf plotrange[1] = minimum(X) end
         if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
         if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
         if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

         W = w[:,:,VarIndex_]

         xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
         yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

         spline = Spline2D(X, Y, W)
         Xi = [i for i in xi, j in yi]
         Yi = [j for i in xi, j in yi]
         wi = spline(Xi[:], Yi[:])
         wi = reshape(wi, size(Xi))'
      end
   end
   return xi, yi, wi
end

# that was pretty easy!
#plot(data[1], head[1], "p")
