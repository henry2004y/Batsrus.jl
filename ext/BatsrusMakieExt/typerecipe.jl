# Type conversion from Batsrus to Makie

"""
Conversion for 1D plots
"""
function Makie.convert_arguments(P::Makie.PointBased, bd::BatsrusIDL, var::String)
   var_ = findindex(bd, var)
   x = bd.x
   y = bd.w[:, var_]

   if hasunit(bd)
      unitx = getunit(bd, bd.head.wname[1])
      unitw = getunit(bd, var)
      if unitx isa UnitfulBatsrus.Unitlike
         x = x .* unitx
      end
      if unitw isa UnitfulBatsrus.Unitlike
         y = y .* unitw
      end
   end

   ([Makie.Point2f(i, j) for (i, j) in zip(x, y)],)
end

"""
Conversion for 2D plots.
"""
function Makie.convert_arguments(P::Makie.GridBased, bd::BatsrusIDL, var::String;
      plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1)
   x, y, w = interp2d(bd, var, plotrange, plotinterval)

   unitx = getunit(bd, bd.head.wname[1])
   unity = getunit(bd, bd.head.wname[2])
   unitw = getunit(bd, var)

   if unitx isa UnitfulBatsrus.Unitlike
      x *= unitx
   end
   if unity isa UnitfulBatsrus.Unitlike
      y *= unity
   end
   if unitw isa UnitfulBatsrus.Unitlike
      w *= unitw
   end

   (x, y, w')
end
