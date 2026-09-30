using LinearAlgebra: tr

# Data manipulation.

"""
    cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
@inline function _check_rectilinear(x, ndim::Int)
    for d in 1:ndim
        DimensionalData.lookup(x, d) isa DimensionalData.NoLookup &&
            error("Selectors are only supported for rectilinear grids (gencoord=false)!")
    end
    return
end

# Convert [lo1 hi1 lo2 hi2 ...] into Between selectors, replacing infinite
# bounds with the endpoints of the corresponding dimension of `data`.
@inline function _resolve_limits(limits, data, ::Val{N}) where {N}
    return ntuple(N) do d
        lo = limits[2d - 1]
        hi = limits[2d]
        lo = isinf(lo) ? first(dims(data, d)) : lo
        hi = isinf(hi) ? last(dims(data, d)) : hi
        return Between(lo, hi)
    end
end

function cutdata(
        bd::BatsrusIDL, var::AbstractString;
        plotrange = [-Inf, Inf, -Inf, Inf], dir::String = "x", sequence::Int = 1
    )
    var_ = findindex(bd, var)

    if dir == "x"
        dim, d1, d2 = 1, 2, 3
    elseif dir == "y"
        dim, d1, d2 = 2, 1, 3
    else
        dim, d1, d2 = 3, 1, 2
    end

    cut1 = selectdim(view(bd.x, :, :, :, d1), dim, sequence)
    cut2 = selectdim(view(bd.x, :, :, :, d2), dim, sequence)
    W = selectdim(view(bd.w, :, :, :, var_), dim, sequence)

    if !all(isinf, plotrange)
        cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
    end

    return cut1, cut2, W
end

@inline function checkvalidlimits(limits, dim::Int = 2)
    return if dim == 2
        if length(limits) != 4
            throw(ArgumentError("Reduction range $limits should be [xmin xmax ymin ymax]!"))
        end

        if limits[1] > limits[2] || limits[3] > limits[4]
            throw(DomainError(limits, "Invalid reduction range!"))
        end
    elseif dim == 3
        if length(limits) != 6
            throw(
                ArgumentError(
                    "Reduction range $limits should be [xmin xmax ymin ymax zmin max]!"
                ),
            )
        end

        if limits[1] > limits[2] || limits[3] > limits[4] || limits[5] > limits[6]
            throw(DomainError(limits, "Invalid reduction range!"))
        end
    end
end

"""
    subsurface(x, y, data, limits)
    subsurface(x, y, u, v, limits)

Extract subset of 2D surface dataset in ndgrid format. See also: [`subvolume`](@ref).
"""
function subsurface(x, y, data, limits)
    checkvalidlimits(limits)
    _check_rectilinear(x, 2)
    selectors = _resolve_limits(limits, data, Val(2))

    # This assumes x and y are the dimensions of data, which they should be for cutdata results
    subdata = data[selectors...]
    # In cutdata, cut1 (x) and cut2 (y) are DimArrays.
    subx = x[selectors...]
    suby = y[selectors...]

    return subx, suby, subdata
end

function subsurface(x, y, u, v, limits)
    checkvalidlimits(limits)
    _check_rectilinear(x, 2)
    selectors = _resolve_limits(limits, u, Val(2))

    return x[selectors...], y[selectors...], u[selectors...], v[selectors...]
end

"""
    subvolume(x, y, z, data, limits)
    subvolume(x, y, z, u, v, w, limits)

Extract subset of 3D dataset in ndgrid format. See also: [`subsurface`](@ref).
"""
function subvolume(x, y, z, data, limits)
    checkvalidlimits(limits, 3)
    _check_rectilinear(x, 3)
    selectors = _resolve_limits(limits, data, Val(3))

    return x[selectors...], y[selectors...], z[selectors...], data[selectors...]
end

function subvolume(x, y, z, u, v, w, limits)
    checkvalidlimits(limits, 3)
    _check_rectilinear(x, 3)
    selectors = _resolve_limits(limits, u, Val(3))

    return x[selectors...], y[selectors...], z[selectors...],
        u[selectors...], v[selectors...], w[selectors...]
end

"""
    getvar(bd::BATS, var::AbstractString) -> Array

Return variable data from string `var`. This is also supported via direct indexing.
Note that the query variable `var` must be in lowercase!

For derived/computed quantities, you can also pass a `Symbol` for a fully
type-stable result:

  - `:b`            — magnetic field magnitude
  - `:b2`           — magnetic field magnitude squared
  - `:e`            — electric field magnitude
  - `:u`            — bulk velocity magnitude
  - `:anisotropy0`  — pressure anisotropy (2D only, species 0)
  - `:anisotropy1`  — pressure anisotropy (2D only, species 1)

# Examples

```julia
bd["rho"]        # direct file variable (string)
bd[:b]           # derived magnitude (symbol, type-stable)
```
"""
function getvar(
        bd::BatsrusIDL{ndim, TV}, var::AbstractString
    ) where {ndim, TV}
    varIndex_ = findindex(bd, var)
    return selectdim(bd.w, ndims(bd.w), varIndex_)
end

"""
Type-stable getvar dispatch via `Val`. The compiler specialises on the symbol
and returns a concretely typed array with no runtime branching inside the loop.
"""
@inline getvar(bd::BatsrusIDL, var::Symbol) = _getvar(bd, Val(var))

# Fallback: treat the symbol as a lowercase string variable name
@inline function _getvar(
        bd::BatsrusIDL{ndim, TV}, ::Val{V}
    ) where {ndim, TV, V}
    varIndex_ = findindex(bd, string(V))
    return selectdim(bd.w, ndims(bd.w), varIndex_)
end

"""
    get_vectors_indices(bd::BatsrusIDL, var::Symbol)

Return indices of vector components for `var`. Supported symbols are `:B`, `:U`, `:E`,
`:U0`, and `:U1`.
"""
@inline get_vectors_indices(bd::BatsrusIDL, var::Symbol) = get_vectors_indices(bd, Val(var))

@inline function get_vectors_indices(bd::BatsrusIDL, ::Val{:B})
    idx = findindex(bd, "bx")
    return (idx, idx + 1, idx + 2)
end

@inline function get_vectors_indices(bd::BatsrusIDL, ::Val{:U})
    idx = findindex(bd, "ux")
    return (idx, idx + 1, idx + 2)
end

@inline function get_vectors_indices(bd::BatsrusIDL, ::Val{:E})
    idx = findindex(bd, "ex")
    return (idx, idx + 1, idx + 2)
end

@inline function get_vectors_indices(bd::BatsrusIDL, ::Val{:U0})
    idx = findindex(bd, "uxs0")
    return (idx, idx + 1, idx + 2)
end

@inline function get_vectors_indices(bd::BatsrusIDL, ::Val{:U1})
    idx = findindex(bd, "uxs1")
    return (idx, idx + 1, idx + 2)
end

@inline function get_vectors_indices(bd::BatsrusIDL, ::Val{:J})
    idx = findindex(bd, "jx")
    return (idx, idx + 1, idx + 2)
end

function get_vectors_indices(bd::BatsrusIDL, ::Val{V}) where {V}
    error("Unknown vector variable $V")
end

"""
    get_vectors(bd::BatsrusIDL, var::Symbol)

Return vector components for `var` as a tuple of arrays.
"""
@inline get_vectors(bd::BatsrusIDL, var::Symbol) = get_vectors(bd, Val(var))

@inline function get_vectors(bd::BatsrusIDL, ::Val{V}) where {V}
    indices = get_vectors_indices(bd, Val(V))
    w = parent(bd.w)
    d = ndims(w)
    return ntuple(i -> selectdim(w, d, indices[i]), length(indices))
end

@inline Base.@propagate_inbounds Base.getindex(bd::BatsrusIDL, var) =
    getvar(bd, var)

"""
    get_timeseries(files::AbstractArray, loc; tstep = 1.0)

Extract plasma moments and EM field from PIC output `files` at `loc` with nearest neighbor.
Currently only works for 2D outputs. If a single point variable is needed, see [`interp1d`](@ref).
"""
# Variables extracted by get_timeseries, in output row order.
const _TIMESERIES_VARS = (
    :rhos0, :rhos1, :uxs0, :uys0, :uzs0, :uxs1, :uys1, :uzs1,
    :pxxs0, :pyys0, :pzzs0, :pxxs1, :pyys1, :pzzs1,
    :bx, :by, :bz, :ex, :ey, :ez,
)

function get_timeseries(files::AbstractArray, loc; tstep = 1.0)
    nfiles = length(files)
    bd = files[1] |> Batsrus.load
    xrange, yrange = get_range(bd)
    trange = range(bd.head.time, step = tstep, length = nfiles)
    @assert xrange[1] ≤ loc[1] ≤ xrange[end] "x location out of range!"
    @assert yrange[1] ≤ loc[2] ≤ yrange[end] "y location out of range!"
    x_ = searchsortedfirst(xrange, loc[1])
    y_ = searchsortedfirst(yrange, loc[2])
    v = zeros(Float32, length(_TIMESERIES_VARS), nfiles)

    @showprogress dt = 1 desc = "Extracting..." for it in eachindex(files)
        bd = files[it] |> Batsrus.load
        for (i, var) in pairs(_TIMESERIES_VARS)
            v[i, it] = bd[var][x_, y_]
        end
    end

    return trange, v
end

function get_timeseries(files::AbstractArray, loc::Vector; tstep = 1.0)
    return get_timeseries(files, SVector{length(loc)}(loc); tstep)
end
