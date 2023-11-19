# Using user recipes from Plots.

using RecipesBase

# Build a recipe which acts on a custom type.
@recipe function f(data::BATLData, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)

   ndim = data.head.ndim

   if startswith(data.head.headline, "normalized")
      hasunits = false
   else
      hasunits = true
      unitw = getunit(data, var)
   end

   if ndim == 1
      VarIndex_ = findindex(data, var)
      if hasunits
         unitx = getunit(data, data.head.variables[1])
         x = data.x .* unitx
         w = data.w
         y = w[:,VarIndex_] .* unitw
      else
         x, w = data.x, data.w
         y = w[:,VarIndex_]
      end

      @series begin
         seriestype --> :path
         x, y
      end
   elseif ndim == 2
      x, y, w = getdata(data, var, plotrange, plotinterval)

      unitx = getunit(data, data.head.variables[1])
      unity = getunit(data, data.head.variables[2])

      if unitx isa UnitfulBatsrus.Unitlike
         x *= unitx
         y *= unity
         w *= unitw
      end

      @series begin
         seriestype --> :contourf  # use := if you want to force it
         x, y, w
      end
   end
end