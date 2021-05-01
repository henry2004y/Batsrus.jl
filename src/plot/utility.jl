# Utility functions for plotting.

"Prepare 2D data arrays for passing to plotting functions."
function getdata(data::Data, var::AbstractString, plotrange, plotinterval, griddim=1)
   @assert data.head.ndim == 2 "data must be in 2D!"

   x, w = data.x, data.w
   ndim = data.head.ndim
   VarIndex_ = findindex(data, var)

   if data.head.gencoord # Generalized coordinates
      X = @view x[:,:,1]
      Y = @view x[:,:,2]
      W = @view w[:,:,VarIndex_]

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
      triang = @views matplotlib.tri.Triangulation(X[:], Y[:])
      interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
      Xi, Yi = meshgrid(xi, yi)
      Wi = interpolator(Xi, Yi)
   else # Cartesian coordinates
      xrange = range(x[1,1,1], x[end,1,1], length=size(x,1))
      yrange = range(x[1,1,2], x[1,end,2], length=size(x,2))
      if all(isinf.(plotrange))
         xi, yi = xrange, yrange
         Xi, Yi = meshgrid(xi, yi)
         Wi = w[:,:,VarIndex_]'
      else
         if plotrange[1] == -Inf plotrange[1] = minimum(xrange) end
         if plotrange[2] ==  Inf plotrange[2] = maximum(xrange) end
         if plotrange[3] == -Inf plotrange[3] = minimum(yrange) end
         if plotrange[4] ==  Inf plotrange[4] = maximum(yrange) end

         xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
         yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

         spline = @views Spline2D(xrange, yrange, w[:,:,VarIndex_])
         Xi, Yi = meshgrid(xi, yi)
         wi = @views spline(Xi[:], Yi[:])
         Wi = reshape(wi, size(Xi))
      end
   end
   if griddim == 1
      return xi, yi, Wi
   else
      return Xi, Yi, Wi
   end
end

"Find variable index in data."
function findindex(data::Data, var::AbstractString)
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
   isnothing(VarIndex_) && error("$(var) not found in file header variables!")
   VarIndex_
end

"Generating consistent 2D arrays for passing to plotting functions."
function meshgrid(x, y)
   X = [x for _ in y, x in x]
   Y = [y for y in y, _ in x]
   X, Y
end