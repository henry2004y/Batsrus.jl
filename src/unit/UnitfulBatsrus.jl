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
@unit AU "AU" AstronomicalUnit (149597870700)*Unitful.m false

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
   "AU" => AU,
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

function _get_typeunit(headline::String)
   if occursin("PIC", headline)
      return "PIC"
   elseif occursin(" AU ", headline)
      return "OUTERHELIO"
   elseif occursin(r"(kg/m3)|(m/s)", headline)
      return "SI"
   elseif occursin(r"(nPa)|( nT )", headline) || headline == "PLANETARY"
      # Some files might have "PLANETARY" in headline but also units?
      # But setunits logic prioritized this.
      return "PLANETARY"
   elseif occursin(r"(dyne)|( G)", headline)
      return "SOLAR"
   else
      return "NORMALIZED"
   end
end

function _get_system_units(typeunit::String)
   # Returns (distance, time, density, velocity, pressure, magnetic field, current density)
   if typeunit == "SI"
      return (Unitful.m, Unitful.s, Unitful.kg / Unitful.m^3, Unitful.m / Unitful.s,
         Unitful.Pa, Unitful.T, Unitful.A / Unitful.m^2)
   elseif typeunit == "CGS"
      # Assuming Unitful.CGS or similar are available/parsed. 
      # Using uparse for safety if symbols aren't exported.
      return (Unitful.cm, Unitful.s, Unitful.g / Unitful.cm^3,
         Unitful.cm / Unitful.s, Unitful.dyn / Unitful.cm^2,
         Unitful.Gauss, Unitful.A / Unitful.m^2) # J in CGS is tricky, defaulting to SI-like or just A/m^2 for validity
   elseif typeunit == "PLANETARY"
      return (UnitfulBatsrus.Re, Unitful.s, UnitfulBatsrus.amucc,
         Unitful.km / Unitful.s, Unitful.nPa, Unitful.nT, Unitful.μA / Unitful.m^2)
   elseif typeunit == "OUTERHELIO"
      return (UnitfulBatsrus.AU, Unitful.s, UnitfulBatsrus.amucc, Unitful.km / Unitful.s,
         Unitful.dyn / Unitful.cm^2, Unitful.nT, Unitful.μA / Unitful.m^2)
   elseif typeunit == "SOLAR"
      return (
         UnitfulBatsrus.RSun, Unitful.s, Unitful.g / Unitful.cm^3, Unitful.km / Unitful.s,
         Unitful.dyn / Unitful.cm^2, Unitful.Gauss, Unitful.μA / Unitful.m^2)
   else # NORMALIZED or PIC
      return (Unitful.NoUnits, Unitful.NoUnits, Unitful.NoUnits, Unitful.NoUnits,
         Unitful.NoUnits, Unitful.NoUnits, Unitful.NoUnits)
   end
end

function _infer_unit_from_name(var::String, sys_units)
   x_u, t_u, rho_u, u_u, p_u, b_u, j_u = sys_units
   v = lowercase(var)
   if startswith(v, "rho") || v == "n"
      return rho_u
   elseif startswith(v, "x") || startswith(v, "y") || startswith(v, "z") ||
          startswith(v, "r")
      return x_u
   elseif startswith(v, "u") || startswith(v, "v")
      return u_u
   elseif startswith(v, "b")
      return b_u
   elseif startswith(v, "p")
      return p_u
   elseif startswith(v, "j")
      return j_u
   else
      return nothing
   end
end

function getunit(bd, var)
   # Batsrus has a bug in the 2D cuts of 3D runs: it always outputs the 3
   # coordinate units in the headline. To work around it, here the index is shifted by 1.
   var_ = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname) + bd.head.ndim + 1
   isnothing(var_) && error("$(var) not found in file header variables!")

   if bd.head.headline in ("normalized variables", "PLANETARY")
      # Known special cases where header string parsing might fail or be implicit
      typeunit = bd.head.headline == "PLANETARY" ? "PLANETARY" :
                 _get_typeunit(bd.head.headline)
      sys_units = _get_system_units(typeunit)
      var_unit = _infer_unit_from_name(var, sys_units)
   else
      var_unit_strs = split(bd.head.headline)
      # Ensure index is within bounds, otherwise fallback to inference
      if var_ <= length(var_unit_strs)
         var_unit = _parse_unit(var_unit_strs[var_])
      else
         # Fallback
         typeunit = _get_typeunit(bd.head.headline)
         sys_units = _get_system_units(typeunit)
         var_unit = _infer_unit_from_name(var, sys_units)
      end
   end

   var_unit
end

function getunits(bd)
   var_unit_strs = split(bd.head.headline)

   # If parsing strategies based on strings is risky, we might want to check for known types first
   typeunit = _get_typeunit(bd.head.headline)
   sys_units = _get_system_units(typeunit)

   if length(var_unit_strs) < length(bd.head.wname)
      # Use inference
      return [_infer_unit_from_name(w, sys_units) for w in bd.head.wname]
   else
      return [_parse_unit(str) for str in var_unit_strs]
   end
end

end
