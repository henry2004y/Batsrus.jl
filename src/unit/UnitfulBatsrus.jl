module UnitfulBatsrus

using Unitful: Unitful
import Unitful: q, c, μ0, ϵ0, k, me, mp
using Unitful: @unit, Unitlike
export getunit, getunits

# lengths
@unit Re "Re" EarthRadii (6378)*Unitful.km false
@unit Rg "Rg" GanymedeRadii (2634)*Unitful.km false # Ganymede's radius
@unit RMercury "RMercury" MercuryRadii (2444)*Unitful.km false # Mercury's radius 
@unit RSun "RSun" SolarRadii (695700)*Unitful.km false # Solar radius

# masses
@unit me_bu "mₑ" ElectronMass me false
@unit mp_bu "mₚ" ProtonMass mp false
@unit u_bu "u" AMU 1*Unitful.u false

# charges
@unit q_bu "q" ElementaryCharge q false
@unit k_bu "k" BoltzmannConstant k false

# densities
@unit amucc "amu/cc" AtomicDensity 1 * Unitful.u/Unitful.cm^3 false

# velocities
@unit kms "km/s" KilometerPerSecond 1 * Unitful.km/Unitful.s false

# current densities
@unit ampm2 "μA/m²" CurrentDensity 1 * Unitful.μA/Unitful.m^2 false

# Others
@unit tm2 "nT/m²" MagneticFieldDivergence 1 * Unitful.nT/Unitful.m^2 false
@unit vm2 "V/m²" ElectricFieldDivergence 1 * Unitful.V/Unitful.m^2 false

__init__() = Unitful.register(UnitfulBatsrus)

function getunit(bd, var)
   # Batsrus has a bug in the 2D cuts of 3D runs: it always outputs the 3
   # coordinate units in the headline. To work around it, here the index is shifted by 1.
   var_ = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname) + bd.head.ndim + 1
   isnothing(var_) && error("$(var) not found in file header variables!")
   if bd.head.headline in ("normalized variables", "PLANETARY")
      var_unit = nothing
   else
      var_unit_strs = split(bd.head.headline)

      if var_unit_strs[var_] == "R"
         var_unit = UnitfulBatsrus.Re
      elseif var_unit_strs[var_] == "Mp/cc"
         var_unit = UnitfulBatsrus.amucc
      elseif var_unit_strs[var_] == "uA/m2"
         var_unit = UnitfulBatsrus.ampm2
      elseif var_unit_strs[var_] == "V/m2"
         var_unit = UnitfulBatsrus.vm2
      else
         try
            var_unit = Unitful.uparse(var_unit_strs[var_])
         catch
            # Fallback for unknown units
            var_unit = nothing
         end
      end
   end

   var_unit
end

function getunits(bd)
   var_unit_strs = split(bd.head.headline)
   var_units = [] # needs to be improved!
   for var_unit_str in var_unit_strs
      if var_unit_str == "R"
         var_unit = UnitfulBatsrus.Re
      elseif var_unit_str == "Mp/cc"
         var_unit = UnitfulBatsrus.amucc
      elseif var_unit_str == "uA/m2"
         var_unit = UnitfulBatsrus.ampm2
      elseif var_unit_str == "V/m2"
         var_unit = UnitfulBatsrus.vm2
      else
         try
            var_unit = Unitful.uparse(var_unit_str)
         catch
            var_unit = nothing
         end
      end
      push!(var_units, var_unit)
   end

   var_units
end

end
