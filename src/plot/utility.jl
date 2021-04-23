# Utility functions for plotting.

"Prepare data for passing to plotting functions."
function getdata(data::Data, var::AbstractString, plotrange, plotinterval)
   x, w = data.x, data.w
   ndim = data.head.ndim
   VarIndex_ = findindex(data, var)

   if data.head.gencoord # Generalized coordinates
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
      Xi = [x for x in xi, _ in yi]
      Yi = [y for _ in xi, y in yi]
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
         Xi = [x for x in xi, _ in yi]
         Yi = [y for _ in xi, y in yi]
         wi = spline(Xi[:], Yi[:])
         wi = reshape(wi, size(Xi))'
      end
   end
   xi, yi, wi
end


"Find variable index in data."
function findindex(data::Data, var::AbstractString)
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
   isnothing(VarIndex_) && error("$(var) not found in file header variables!")
   VarIndex_
end
