# Utility functions for plotting.

"""
    getdata2d(bd::BATLData, var::AbstractString, plotrange=[-Inf, Inf, -Inf, Inf],
       plotinterval=Inf; innermask=false)

Return 2D slices of data `var` from `bd`. If `plotrange` is not set, output data resolution
is the same as the original. If `innermask==true`, then the inner boundary cells are set to
NaN.
"""
function getdata2d(bd::BATLData, var::AbstractString,
   plotrange::Vector=[-Inf, Inf, -Inf, Inf], plotinterval::Real=Inf; innermask::Bool=false)
   x, w, ndim = bd.x, bd.w, bd.head.ndim
   @assert ndim == 2 "data must be in 2D!"

   varIndex_ = findindex(bd, var)

   if bd.head.gencoord # Generalized coordinates
      X, Y = eachslice(x, dims=3)
      W = @view w[:,:,varIndex_]

      adjust_plotrange!(plotrange, extrema(X), extrema(Y))
      # Set a heuristic value if not set
      if isinf(plotinterval)
         plotinterval = (plotrange[2] - plotrange[1]) / size(X, 1)
      end
      xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
      yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
      # Perform linear interpolation of the data (x,y) on grid(xi,yi)
      triang = @views matplotlib.tri.Triangulation(X[:], Y[:])
      interpolator = @views matplotlib.tri.LinearTriInterpolator(triang, W[:])
      Xi, Yi = meshgrid(xi, yi)
      Wi = interpolator(Xi, Yi)
   else # Cartesian coordinates
      xrange = range(x[1,1,1], x[end,1,1], length=size(x,1))
      yrange = range(x[1,1,2], x[1,end,2], length=size(x,2))
      if all(isinf.(plotrange))
         xi, yi = xrange, yrange
         Wi = w[:,:,varIndex_]' # Matplotlib does not accept view!
      else
         adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

         if isinf(plotinterval)
            xi = range(plotrange[1], stop=plotrange[2], step=xrange[2] - xrange[1])
            yi = range(plotrange[3], stop=plotrange[4], step=yrange[2] - yrange[1])
         else
            xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
            yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
         end

         interp = @views cubic_spline_interpolation((xrange, yrange), w[:,:,varIndex_])
         Wi = [interp(i, j) for j in yi, i in xi]
      end
   end

   # Mask a circle at the inner boundary
   if innermask
      varIndex_ = findlast(x->x=="rbody", bd.head.variables)
      isnothing(varIndex_) && error("rbody not found in file header parameters!")
      ParamIndex_ = varIndex_ - ndim - bd.head.nw
      @inbounds for i in CartesianIndices(Wi)
         if xi[i[1]]^2 + yi[i[2]]^2 < bd.head.eqpar[ParamIndex_]^2
            Wi[i] = NaN
         end
      end
   end

   xi, yi, Wi
end

"Find variable index in the BATSRUS data."
function findindex(bd::BATLData, var::AbstractString)
   varIndex_ = findfirst(x->x==lowercase(var), lowercase.(bd.head.wnames))
   isnothing(varIndex_) && error("$(var) not found in file header variables!")

   varIndex_
end

"Generating consistent 2D arrays for passing to plotting functions."
function meshgrid(x, y)
   X = [x for _ in y, x in x]
   Y = [y for y in y, _ in x]

   X, Y
end

@inline function hasunit(bd::BATLData)
   if startswith(bd.head.headline, "normalized")
      return false
   else
      return true
   end
end

"Adjust 2D plot ranges."
function adjust_plotrange!(plotrange, xlimit, ylimit)
   plotrange[1] = ifelse(isinf(plotrange[1]), xlimit[1], plotrange[1])
   plotrange[2] = ifelse(isinf(plotrange[2]), xlimit[2], plotrange[2])
   plotrange[3] = ifelse(isinf(plotrange[3]), ylimit[1], plotrange[3])
   plotrange[4] = ifelse(isinf(plotrange[4]), ylimit[2], plotrange[4])

   return
end