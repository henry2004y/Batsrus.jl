using LinearAlgebra: tr

# Data manipulation.

"""
    cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(
        bd::BatsrusIDL, var::AbstractString;
        plotrange = [-Inf, Inf, -Inf, Inf], dir::String = "x", sequence::Int = 1
    )
    var_ = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname)
    isempty(var_) && error("$(var) not found in header variables!")

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

    if !all(isinf.(plotrange))
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

    if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
            DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup
        error("Selectors are only supported for rectilinear grids (gencoord=false)!")
    end

    xmin, xmax, ymin, ymax = limits
    xmin = isinf(xmin) ? first(dims(data, 1)) : xmin
    xmax = isinf(xmax) ? last(dims(data, 1)) : xmax
    ymin = isinf(ymin) ? first(dims(data, 2)) : ymin
    ymax = isinf(ymax) ? last(dims(data, 2)) : ymax

    selectors = (Between(xmin, xmax), Between(ymin, ymax))
    # This assumes x and y are the dimensions of data, which they should be for cutdata results
    subdata = data[selectors...]
    # subx and suby should also be sliced if they are DimArrays, or we can just return the dims from subdata
    # In cutdata, cut1 (x) and cut2 (y) are DimArrays.
    subx = x[selectors...]
    suby = y[selectors...]

    return subx, suby, subdata
end

function subsurface(x, y, u, v, limits)
    checkvalidlimits(limits)

    if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
            DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup
        error("Selectors are only supported for rectilinear grids (gencoord=false)!")
    end

    xmin, xmax, ymin, ymax = limits
    xmin = isinf(xmin) ? first(dims(u, 1)) : xmin
    xmax = isinf(xmax) ? last(dims(u, 1)) : xmax
    ymin = isinf(ymin) ? first(dims(u, 2)) : ymin
    ymax = isinf(ymax) ? last(dims(u, 2)) : ymax

    selectors = (Between(xmin, xmax), Between(ymin, ymax))
    newu = u[selectors...]
    newv = v[selectors...]

    subx = x[selectors...]
    suby = y[selectors...]

    return subx, suby, newu, newv
end

"""
    subvolume(x, y, z, data, limits)
    subvolume(x, y, z, u, v, w, limits)

Extract subset of 3D dataset in ndgrid format. See also: [`subsurface`](@ref).
"""
function subvolume(x, y, z, data, limits)
    checkvalidlimits(limits, 3)

    if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
            DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup ||
            DimensionalData.lookup(x, 3) isa DimensionalData.NoLookup
        error("Selectors are only supported for rectilinear grids (gencoord=false)!")
    end

    xmin, xmax, ymin, ymax, zmin, zmax = limits
    xmin = isinf(xmin) ? first(dims(data, 1)) : xmin
    xmax = isinf(xmax) ? last(dims(data, 1)) : xmax
    ymin = isinf(ymin) ? first(dims(data, 2)) : ymin
    ymax = isinf(ymax) ? last(dims(data, 2)) : ymax
    zmin = isinf(zmin) ? first(dims(data, 3)) : zmin
    zmax = isinf(zmax) ? last(dims(data, 3)) : zmax

    selectors = (Between(xmin, xmax), Between(ymin, ymax), Between(zmin, zmax))
    subdata = data[selectors...]
    subx = x[selectors...]
    suby = y[selectors...]
    subz = z[selectors...]

    return subx, suby, subz, subdata
end

function subvolume(x, y, z, u, v, w, limits)
    checkvalidlimits(limits, 3)

    if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
            DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup ||
            DimensionalData.lookup(x, 3) isa DimensionalData.NoLookup
        error("Selectors are only supported for rectilinear grids (gencoord=false)!")
    end

    xmin, xmax, ymin, ymax, zmin, zmax = limits
    xmin = isinf(xmin) ? first(dims(u, 1)) : xmin
    xmax = isinf(xmax) ? last(dims(u, 1)) : xmax
    ymin = isinf(ymin) ? first(dims(u, 2)) : ymin
    ymax = isinf(ymax) ? last(dims(u, 2)) : ymax
    zmin = isinf(zmin) ? first(dims(u, 3)) : zmin
    zmax = isinf(zmax) ? last(dims(u, 3)) : zmax

    selectors = (Between(xmin, xmax), Between(ymin, ymax), Between(zmin, zmax))
    newu = u[selectors...]
    newv = v[selectors...]
    neww = w[selectors...]

    subx = x[selectors...]
    suby = y[selectors...]
    subz = z[selectors...]

    return subx, suby, subz, newu, newv, neww
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
function get_vectors_indices(bd::BatsrusIDL, var::Symbol)
    vnames = if var == :B
        ["bx", "by", "bz"]
    elseif var == :U
        ["ux", "uy", "uz"]
    elseif var == :E
        ["ex", "ey", "ez"]
    elseif var == :U0
        ["uxs0", "uys0", "uzs0"]
    elseif var == :U1
        ["uxs1", "uys1", "uzs1"]
    else
        error("Unknown vector variable $var")
    end
    return ntuple(i -> findindex(bd, vnames[i]), length(vnames))
end

"""
    get_vectors(bd::BatsrusIDL, var::Symbol)

Return vector components for `var` as a tuple of arrays.
"""
function get_vectors(bd::BatsrusIDL, var::Symbol)
    indices = get_vectors_indices(bd, var)
    return ntuple(i -> selectdim(bd.w, ndims(bd.w), indices[i]), length(indices))
end

@inline Base.@propagate_inbounds Base.getindex(bd::BatsrusIDL, var) =
    getvar(bd, var)

"""
    get_timeseries(files::AbstractArray, loc; tstep = 1.0)

Extract plasma moments and EM field from PIC output `files` at `loc` with nearest neighbor.
Currently only works for 2D outputs. If a single point variable is needed, see [`interp1d`](@ref).
"""
function get_timeseries(files::AbstractArray, loc; tstep = 1.0)
    nfiles = length(files)
    bd = files[1] |> Batsrus.load
    xrange, yrange = get_range(bd)
    trange = range(bd.head.time, step = tstep, length = nfiles)
    @assert xrange[1] ≤ loc[1] ≤ xrange[end] "x location out of range!"
    @assert yrange[1] ≤ loc[2] ≤ yrange[end] "y location out of range!"
    x_ = searchsortedfirst(xrange, loc[1])
    y_ = searchsortedfirst(yrange, loc[2])
    v = zeros(Float32, 20, nfiles)

    @showprogress dt = 1 desc = "Extracting..." for it in eachindex(files)
        bd = files[it] |> Batsrus.load
        v[1, it] = bd[:rhos0][x_, y_]
        v[2, it] = bd[:rhos1][x_, y_]
        v[3, it] = bd[:uxs0][x_, y_]
        v[4, it] = bd[:uys0][x_, y_]
        v[5, it] = bd[:uzs0][x_, y_]
        v[6, it] = bd[:uxs1][x_, y_]
        v[7, it] = bd[:uys1][x_, y_]
        v[8, it] = bd[:uzs1][x_, y_]
        v[9, it] = bd[:pxxs0][x_, y_]
        v[10, it] = bd[:pyys0][x_, y_]
        v[11, it] = bd[:pzzs0][x_, y_]
        v[12, it] = bd[:pxxs1][x_, y_]
        v[13, it] = bd[:pyys1][x_, y_]
        v[14, it] = bd[:pzzs1][x_, y_]
        v[15, it] = bd[:bx][x_, y_]
        v[16, it] = bd[:by][x_, y_]
        v[17, it] = bd[:bz][x_, y_]
        v[18, it] = bd[:ex][x_, y_]
        v[19, it] = bd[:ey][x_, y_]
        v[20, it] = bd[:ez][x_, y_]
    end

    return trange, v
end

function get_timeseries(files::AbstractArray, loc::Vector; tstep = 1.0)
    return get_timeseries(files, SVector{length(loc)}(loc); tstep)
end
