module FlexUnitsBatsrus

using FlexUnits
using FlexUnits.UnitRegistry
export getunit, getunits, me, mp, q, k, μ0, ε0, R, Ry

const _me = Ref{Any}()
const _mp = Ref{Any}()
const _q = Ref{Any}()
const _k = Ref{Any}()
const _μ0 = Ref{Any}()
const _ε0 = Ref{Any}()
const _R = Ref{Any}()
const _Ry = Ref{Any}()

# Accessors for constants to be used as values
me() = _me[]
mp() = _mp[]
q() = _q[]
k() = _k[]
μ0() = _μ0[]
ε0() = _ε0[]
R() = _R[]
Ry() = _Ry[]

const _UNIT_MAP = Dict{String, Any}()

function __init__()
    # Skip registration during precompilation of this module or its extensions
    # because uparse/register_unit use eval, which breaks incremental compilation.
    if Base.generating_output()
        return
    end
    # register custom units
    # Re: Earth Radii
    register_unit("Re" => qparse("6378km"))
    # Rg: Ganymede Radii
    register_unit("Rg" => qparse("2634km"))
    # RMercury: Mercury Radii
    register_unit("RMercury" => qparse("2444km"))
    # RSun: Solar Radii
    register_unit("RSun" => qparse("695700km"))
    # AU: Astronomical Unit
    register_unit("AU" => qparse("149597870700m"))
    # u: unified atomic mass unit
    register_unit("u" => qparse("1.66053906660e-27kg"))
    # nT: nanoTesla
    register_unit("nT" => 1e-9 * qparse("T"))
    # nPa: nanoPascal
    register_unit("nPa" => 1e-9 * qparse("Pa"))
    # eV: electronvolt
    register_unit("eV" => 1.602176634e-19 * qparse("J"))
    # Ry: Rydberg
    register_unit("Ry" => 13.605693122994 * 1.602176634e-19 * qparse("J"))

    # densities
    # amucc: amu/cc
    register_unit("amucc" => qparse("u/cm^3"))

    # velocities
    # kms: km/s
    register_unit("kms" => qparse("km/s"))

    # current densities
    # ampm2: μA/m²
    register_unit("ampm2" => qparse("μA/m^2"))

    # Others
    # tm2: nT/m²
    register_unit("tm2" => qparse("nT/m^2"))
    # vm2: V/m²
    register_unit("vm2" => qparse("V/m^2"))

    # Physical constants
    _me[] = 9.1093837015e-31 * qparse("kg")
    _mp[] = 1.67262192369e-27 * qparse("kg")
    _q[] = 1.602176634e-19 * qparse("C")
    _k[] = 1.380649e-23 * qparse("J/K")
    _μ0[] = 4π * 1e-7 * qparse("T*m/A")
    _ε0[] = 8.8541878128e-12 * qparse("F/m")
    _R[] = 8.31446261815324 * qparse("J/(mol*K)")
    _Ry[] = 13.605693122994 * 1.602176634e-19 * qparse("J")

    # Unit Map
    _UNIT_MAP["R"] = qparse("Re").unit
    _UNIT_MAP["AU"] = qparse("AU").unit
    _UNIT_MAP["Mp/cc"] = qparse("amucc").unit
    _UNIT_MAP["uA/m2"] = qparse("ampm2").unit
    _UNIT_MAP["V/m2"] = qparse("vm2").unit
end

function _parse_unit(unit_str::AbstractString)
    # First, check our custom mapping.
    unit = get(_UNIT_MAP, unit_str, nothing)
    if !isnothing(unit)
        return unit
    end

    # FlexUnits doesn't have a direct equivalent of Unitful.uparse for arbitrary strings
    # but we can try to use its qparse function and extract the unit.
    try
        # Handle some common non-standard strings
        if unit_str == "kg/m3"
            return qparse("kg/m^3").unit
        elseif unit_str == "dyne"
            return qparse("dyn").unit
        end
        return qparse(String(unit_str)).unit
    catch e
        return nothing
    end
end

function _get_typeunit(headline::AbstractString)
    if occursin("PIC", headline)
        return "PIC"
    elseif occursin(" AU ", headline)
        return "OUTERHELIO"
    elseif occursin(r"(kg/m3)|(m/s)", headline)
        return "SI"
    elseif occursin(r"(nPa)|( nT )", headline) || headline == "PLANETARY"
        return "PLANETARY"
    elseif occursin(r"(dyne)|( G)", headline)
        return "SOLAR"
    else
        return "NORMALIZED"
    end
end

function _get_system_units(typeunit::AbstractString)
    if typeunit == "SI"
        return (
            qparse("m").unit, qparse("s").unit, qparse("kg/m^3").unit, qparse("m/s").unit, qparse("Pa").unit, qparse("T").unit, qparse("A/m^2").unit
        )
    elseif typeunit == "CGS"
        return (
            qparse("cm").unit, qparse("s").unit, qparse("g/cm^3").unit, qparse("cm/s").unit, qparse("dyn/cm^2").unit, qparse("G").unit, qparse("A/m^2").unit
        )
    elseif typeunit == "PLANETARY"
        return (
            qparse("Re").unit, qparse("s").unit, qparse("amucc").unit, qparse("km/s").unit, qparse("nPa").unit, qparse("nT").unit, qparse("μA/m^2").unit
        )
    elseif typeunit == "OUTERHELIO"
        return (
            qparse("AU").unit, qparse("s").unit, qparse("amucc").unit, qparse("km/s").unit, qparse("dyn/cm^2").unit, qparse("nT").unit, qparse("μA/m^2").unit
        )
    elseif typeunit == "SOLAR"
        return (
            qparse("RSun").unit, qparse("s").unit, qparse("g/cm^3").unit, qparse("km/s").unit, qparse("dyn/cm^2").unit, qparse("G").unit, qparse("μA/m^2").unit
        )
    else # NORMALIZED or PIC
        return (nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end

function _infer_unit_from_name(var::AbstractString, sys_units)
    x_u, t_u, rho_u, u_u, p_u, b_u, j_u = sys_units
    v = lowercase(var)
    if isnothing(x_u) return nothing end

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
    idx = findfirst(x -> lowercase(x) == lowercase(var), bd.head.coord)
    if isnothing(idx)
        idx = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname)
        isnothing(idx) && error("$(var) not found in file header variables!")
        var_ = idx + bd.head.ndim + 1
    else
        var_ = idx
    end

    if bd.head.headline in ("normalized variables", "PLANETARY")
        typeunit = bd.head.headline == "PLANETARY" ? "PLANETARY" :
            _get_typeunit(bd.head.headline)
        sys_units = _get_system_units(typeunit)
        var_unit = _infer_unit_from_name(var, sys_units)
    else
        var_unit_strs = split(bd.head.headline)
        if var_ <= length(var_unit_strs)
            var_unit = _parse_unit(var_unit_strs[var_])
        else
            typeunit = _get_typeunit(bd.head.headline)
            sys_units = _get_system_units(typeunit)
            var_unit = _infer_unit_from_name(var, sys_units)
        end
    end

    return var_unit
end

function getunits(bd)
    var_unit_strs = split(bd.head.headline)
    typeunit = _get_typeunit(bd.head.headline)
    sys_units = _get_system_units(typeunit)

    if length(var_unit_strs) < length(bd.head.wname)
        return [_infer_unit_from_name(w, sys_units) for w in bd.head.wname]
    else
        return [_parse_unit(str) for str in var_unit_strs]
    end
end

end
