# Using user recipes from Plots.

using RecipesBase

# Build a recipe which acts on a custom type.
@recipe function f(bd::BATLData, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)

   ndim = bd.head.ndim

   hasunits = hasunit(bd)

   if ndim == 1
      VarIndex_ = findindex(bd, var)
      if hasunits
         unitx = getunit(bd, bd.head.variables[1])
         unitw = getunit(bd, var)
         x = bd.x .* unitx
         y = bd.w[:,VarIndex_] .* unitw
      else
         x = bd.x
         y = @view bd.w[:,VarIndex_]
      end

      @series begin
         seriestype --> :path
         x, y
      end
   elseif ndim == 2
      x, y, w = getdata(bd, var, plotrange, plotinterval)

      unitx = getunit(bd, bd.head.variables[1])
      unity = getunit(bd, bd.head.variables[2])
      unitw = getunit(bd, var)

      if unitx isa UnitfulBatsrus.Unitlike
         x *= unitx
         y *= unity
         w *= unitw
      end

      @series begin
         seriestype --> :contourf  # use := if you want to force it
         x, y, w'
      end
   end
end