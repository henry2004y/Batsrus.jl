# Using user recipes from Plots.

using RecipesBase

# Build a recipe which acts on a custom type.
@recipe function f(bd::AbstractBATS, var::AbstractString;
      plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1)
   hasunits = hasunit(bd)
   if bd.head.ndim == 1
      if hasunits
         unitx = getunit(bd, bd.head.wname[1])
         unitw = getunit(bd, var)
         x = bd.x .* unitx
         y = getview(bd, var) .* unitw
      else
         x = bd.x
         y = getview(bd, var)
      end

      @series begin
         seriestype --> :path
         x, y
      end
   elseif bd.head.ndim == 2
      x, y, w = interp2d(bd, var, plotrange, plotinterval)
      if hasunits
         unitx = getunit(bd, bd.head.wname[1])
         unity = getunit(bd, bd.head.wname[2])
         unitw = getunit(bd, var)

         if unitx isa UnitfulBatsrus.Unitlike
            x *= unitx
            y *= unity
            w *= unitw
         end
      end

      @series begin
         seriestype --> :contourf  # use := if you want to force it
         x, y, w'
      end
   end
end
