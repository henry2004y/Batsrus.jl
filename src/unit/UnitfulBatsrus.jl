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

const _UNIT_MAP = Dict(
   "R" => Re,
   "Mp/cc" => amucc,
   "uA/m2" => ampm2,
   "V/m2" => vm2
)

__init__() = Unitful.register(UnitfulBatsrus)

function _parse_unit(unit_str::AbstractString)
   # First, check our custom mapping.
   unit = get(_UNIT_MAP, unit_str, nothing)
   if !isnothing(unit)
      return unit
   end

   # If not found, try standard parsing.
   try
      return Unitful.uparse(unit_str)
   catch
      # Return nothing if parsing fails.
      return nothing
   end
end

function getunit(bd, var)
   # Batsrus has a bug in the 2D cuts of 3D runs: it always outputs the 3
   # coordinate units in the headline. To work around it, here the index is shifted by 1.
   var_ = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname) + bd.head.ndim + 1
   isnothing(var_) && error("$(var) not found in file header variables!")
   if bd.head.headline in ("normalized variables", "PLANETARY")
      var_unit = nothing
   else
      var_unit_strs = split(bd.head.headline)
      var_unit = _parse_unit(var_unit_strs[var_])
   end

   var_unit
end

function getunits(bd)
   var_unit_strs = split(bd.head.headline)
   var_units = [] # needs to be improved!
   for var_unit_str in var_unit_strs
      var_unit = _parse_unit(var_unit_str)
      push!(var_units, var_unit)
   end

   var_units
end

end
