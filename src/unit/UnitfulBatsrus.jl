module UnitfulBatsrus

import Unitful
import Unitful: q, c, μ0, ϵ0, k, me, mp
using Unitful: @unit, Unitlike
export @bu_str
export getunit, getunits


# lengths
@unit R_bu  "Rₑ"     EarthRadii        (6378)*Unitful.km        false
@unit Rg_bu  "Rg"    GanymedeRadii     (2634)*Unitful.km        false # Ganymede's radius
@unit Rm_bu  "Rm"    MercuryRadii      (2444)*Unitful.km        false # Mercury's radius 

# masses
@unit me_bu "mₑ"     ElectronMass      me                       false
@unit mp_bu "mₚ"     ProtonMass        mp                       false
@unit u_bu  "u"      AMU               1*Unitful.u              false

# charges
@unit q_bu  "q"      ElementaryCharge  q                        false
@unit k_bu  "k"      BoltzmannConstant k                        false

# densities
@unit amucc_bu "amu/cc" AtomicDensity  1*Unitful.u/Unitful.cm^3 false

# velocities
@unit V_bu  "V"      Speed             1*Unitful.km/Unitful.s   false

# current densities
@unit ampm2_bu  "μA/m²" CurrentDensity  1*Unitful.μA/Unitful.m^2 false

# Others
@unit tm2_bu "nT/m²" MagneticFieldDivergence 1*Unitful.nT/Unitful.m^2 false
@unit vm2_bu "V/m²"  ElectricFieldDivergence 1*Unitful.V/Unitful.m^2  false

include("batsrusmacro.jl")

# Some gymnastics required here because if we precompile, we cannot add to
# Unitful.basefactors at compile time and expect the changes to persist to runtime.
const localunits = Unitful.basefactors
function __init__()
   merge!(Unitful.basefactors, localunits)
   Unitful.register(UnitfulBatsrus)
end


function getunit(bd, var)
   # Batsrus has a bug in the 2D cuts of 3D runs: it always outputs the 3
   # coordinate units in the headline. To work around it, here the index is shifted by 1.
   var_ = findfirst(x->lowercase(x)==lowercase(var), bd.head.wname) + bd.head.ndim + 1
   isnothing(var_) && error("$(var) not found in file header variables!")
   if bd.head.headline in ("normalized variables", "PLANETARY")
      var_unit = nothing 
   else
      var_unit_strs = split(bd.head.headline)
 
      if var_unit_strs[var_] == "R"
         var_unit = bu"R"
      elseif var_unit_strs[var_] == "Mp/cc"
         var_unit = bu"amucc"
      elseif var_unit_strs[var_] == "uA/m2"
         var_unit = bu"ampm2"
      elseif var_unit_strs[var_] == "V/m2"
         var_unit = bu"vm2"
      else
         var_unit = Unitful.uparse(var_unit_strs[var_])
      end
   end
   
   var_unit
end
 
function getunits(bd)
   var_unit_strs = split(bd.head.headline)
   var_units = [] # needs to be improved!
   for var_unit_str in var_unit_strs
      if var_unit_str == "R"
         var_unit = bu"R"
      elseif var_unit_str == "Mp/cc"
         var_unit = bu"amucc"
      elseif var_unit_str == "uA/m2"
         var_unit = bu"ampm2"
      elseif var_unit_str == "V/m2"
         var_unit = bu"vm2"
      else
         var_unit = UnitfulBatsrus.Unitful.uparse(var_unit_str)
      end
      push!(var_units, var_unit)
   end

   var_units
end

end