# Type conversion from Batsrus to Makie

"Conversion for 1D plots"
function Makie.convert_arguments(P::Makie.PointBased, bd::BATS, var::String)
   var_ = findindex(bd, var)
   if hasunit(bd)
      unitx = getunit(bd, bd.head.variables[1])
      unitw = getunit(bd, var)
      x = bd.x .* unitx
      y = bd.w[:,var_] .* unitw
   else
      x = bd.x
      y = bd.w[:,var_]
   end

   ([Makie.Point2f(i, j) for (i, j) in zip(x, y)],)
end

"Conversion for 2D plots."
function Makie.convert_arguments(P::Makie.GridBased, bd::BATS, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)
   x, y, w = interp2d(bd, var, plotrange, plotinterval)

   unitx = getunit(bd, bd.head.variables[1])
   unity = getunit(bd, bd.head.variables[2])
   unitw = getunit(bd, var)

   if unitx isa UnitfulBatsrus.Unitlike
      x *= unitx
      y *= unity
      w *= unitw
   end

   (x, y, w')
end
