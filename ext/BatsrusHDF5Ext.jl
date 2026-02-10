module BatsrusHDF5Ext

using Batsrus
using HDF5
using StaticArrays: SVector, MVector

import Batsrus: BatsrusHDF5File, extract_var

export BatsrusHDF5Uniform

"""
BATSRUS hdf5 file wrapper.

The data are stored in blocks, i.e., each field component is stored in a 4D array in the
order (iblock, iz, iy, ix). This is a generic wrapper and does not assume grid type, i.e.,
uniform, stretched nonuniform, or AMR, etc. Classes to handle data with different grids can
be derived from this class.
"""
struct HDF5Common{TI <: Signed, TF <: AbstractFloat} <: BatsrusHDF5File
    fid::HDF5.File
    version::TI
    "Type of geometry. 0 for Cartesian."
    geometry::TI
    "Saved snapshot timestamp."
    time::TF
    "Minimum coordinates for the whole grid."
    coordmin::SVector{3, TF}
    "Maximum coordinates for the whole grid."
    coordmax::SVector{3, TF}
    "Saved snapshot simulation timestep."
    timestep::TI
    "True dimension of the data despite the stored data dimension."
    ndim::TI
    "Number of cell in each block along x, y, z."
    ncb::SVector{3, TI}
    "If boundary condition is periodic."
    isperiodic::SVector{3, Bool}
    "Vector of non-singleton dimensions."
    multi_cell_dims::SVector{3, Bool}
    "Lengths along each direction for the whole grid."
    extent::SVector{3, TF}

    function HDF5Common(filename::AbstractString)
        fid = h5open(filename, "r")
        meta_int = read(fid["Integer Plot Metadata"])::Vector{Int32}
        meta_real = read(fid["Real Plot Metadata"])::Vector{<:AbstractFloat}

        version = meta_int[1]
        geometry = meta_int[11]
        time = meta_real[1]
        coordmin = SVector(meta_real[2], meta_real[4], meta_real[6])
        coordmax = SVector(meta_real[3], meta_real[5], meta_real[7])
        timestep = meta_int[2]
        ndim = meta_int[3]
        ncb = SVector(meta_int[8], meta_int[9], meta_int[10])
        isperiodic = SVector(meta_int[12], meta_int[13], meta_int[14])

        multi_cell_dims = ncb .!= 0
        extent = coordmax .- coordmin

        return new{eltype(meta_int), eltype(meta_real)}(
            fid, version, geometry, time, coordmin,
            coordmax, timestep, ndim, ncb, isperiodic, multi_cell_dims, extent
        )
    end
end

"""
BATSRUS HDF5 file with uniform Cartesian mesh.
"""
struct _BatsrusHDF5Uniform{TI, TF} <: BatsrusHDF5File
    common::HDF5Common{TI, TF}
    "Numbers of cells along each direction"
    nc::SVector{3, TI}
    "Number of blocks along each direction"
    nb::SVector{3, TI}
    "Grid resolution along each direction"
    dcoord::SVector{3, TF}
    "Block length along each direction"
    dblock::SVector{3, TF}

    function _BatsrusHDF5Uniform(filename::AbstractString)
        bf = HDF5Common(filename)

        nc = MVector{3, Int32}(1, 1, 1)
        try
            extent::Matrix{Int32} = bf.fid["MinLogicalExtents"] |> read
            nc[bf.multi_cell_dims] = extent[:, end] + bf.ncb[bf.multi_cell_dims]
        catch
            iCoord_DB::Matrix{Int32} = bf.fid["iCoord_DB"] |> read
            nc[bf.multi_cell_dims] = iCoord_DB[:, end] + bf.ncb[bf.multi_cell_dims]
        end
        nb = SVector{3}(nc .÷ bf.ncb)

        dcoord = SVector{3}(bf.extent ./ nc)
        dblock = SVector{3}(bf.extent ./ nb)

        TI, TF = findparam(bf)

        return new{TI, TF}(bf, SVector{3}(nc), nb, dcoord, dblock)
    end
end

Batsrus.BatsrusHDF5Uniform(filename::AbstractString) = _BatsrusHDF5Uniform(filename)

findparam(::HDF5Common{TI, TF}) where {TI, TF} = (TI, TF)

function Base.show(io::IO, file::_BatsrusHDF5Uniform)
    println(io, "Dimension                : ", file.common.ndim)
    println(io, "Mesh coordmin            : ", file.common.coordmin)
    println(io, "Mesh coordmax            : ", file.common.coordmax)
    println(io, "Number of blocks         : ", file.nb)
    println(io, "Number of cells per block: ", file.common.ncb)
    println(io, "Grid resolution          : ", file.dcoord)
    println(io, "Time                     : ", file.common.time)
    vars = HDF5.keys(file.common.fid)
    idBegin_ = findfirst(x -> x == "bounding box", vars) + 1
    idEnd_ = findfirst(x -> endswith(x, "Ext"), vars) - 1
    return println(io, "Variables                : ", vars[idBegin_:idEnd_])
end

"""
    prep_extract(file::BatsrusHDF5Uniform, vmin=-Inf, vmax=Inf, step=1)

Get info for data extraction in 1D.

# Keywords

  - `vmin, vmax`: requested coordinate range (corner values).
  - `step`: stride.

# Returns:

  - `gslc`: global slice.
  - `vmin_new, vmax_new`: adjusted coordinate range (corner values).
  - `ibmin:ibmax`: block range.
"""
function prep_extract(
        file::_BatsrusHDF5Uniform;
        dim::Int = 1, vmin = -Inf32, vmax = Inf32, step::Int = 1
    )
    vmin = isinf(vmin) ? file.common.coordmin[dim] : max(file.common.coordmin[dim], vmin)
    vmax = isinf(vmax) ? file.common.coordmax[dim] : min(file.common.coordmax[dim], vmax)
    # global slice and adjusted coordinate bounds
    gslc, vmin_new, vmax_new = prepslice(file; dim, vmin, vmax, step)
    # blocks involved
    ibmin = gslc.start ÷ file.common.ncb[dim]
    ibmax = gslc.stop ÷ file.common.ncb[dim] - 1
    # If a global slice's stop lies in the interior of a block, this block is included.
    if gslc.stop % file.common.ncb[dim] > 0
        ibmax += 1
    end
    ibmax = max(ibmax, ibmin)

    return gslc, vmin_new, vmax_new, ibmin:ibmax
end

"""
Return lower corner index.
"""
function coord2index(file::_BatsrusHDF5Uniform{TI, TF}, dim::Int, x::Real) where {TI, TF}
    if file.nc[dim] == 1
        return 1
    else
        return ceil(Int, (x - file.common.coordmin[dim]) / file.dcoord[dim] + TF(0.1))
    end
end

"""
    prepslice(file::BatsrusHDF5Uniform; dim::Int, vmin, vmax, step=1)

Return range that covers [`vmin`, `vmax`) along dimension `dim`.

# Returns

  - `slc_new`: trimed slice. If the object's Nx==1, then 1:1 will be returned.
  - `xl_new`: adjusted lower corner coordinate matching `slc_new.start`.
  - `xu_new`: adjusted lower corner coordinate matching `slc_new.stop`.
"""
function prepslice(file::_BatsrusHDF5Uniform; dim::Int = 1, vmin, vmax, step::Int = 1)
    start = coord2index(file, dim, vmin)
    stop = max(coord2index(file, dim, vmax), start)

    slc_new = trimslice(start, stop, step, file.nc[dim])
    xl_new = file.common.coordmin[dim] + (slc_new.start - 1) * file.dcoord[dim]
    xu_new = file.common.coordmin[dim] + slc_new.stop * file.dcoord[dim]

    return slc_new, xl_new, xu_new
end

"""
    trimslice(start, stop, step, stop_max)

Set slice's start to be nonnegative and start/stop to be within bound.
Reverse slicing is not handled.
"""
function trimslice(start, stop, step, stop_max)
    if start < 1
        start += (-start ÷ step) * step
        if start < 1
            start += step
        end
    end
    start = min(start, stop_max)
    stop = min(max(stop, start), stop_max)

    return start:step:stop
end

"""
    prep_extract_per_block(file::BatsrusHDF5Uniform, dim, gslc, ib)

Get info for data extraction on a single block.

# Arguments

  - `gslc`: global slice from prep_extract.
  - `ib::Int`: block index.

# Returns

  - `lslc`: range to be used on the current block.
  - `ix0:ix1`: index range in the global array.
"""
@inline function prep_extract_per_block(
        file::_BatsrusHDF5Uniform, dim::Int,
        gslc::OrdinalRange, ib::Int
    )
    # compute local slice on this block
    lslc = global_slice_to_local_slice(file, dim, gslc, ib)
    # compute the index bounds in the output array
    offset = lslc.start + file.common.ncb[dim] * ib - gslc.start
    ix0 = offset ÷ gslc.step + 1
    ix1 = ix0 + length(lslc) - 1

    return lslc, ix0:ix1
end

"""
    global_slice_to_local_slice(file::BatsrusHDF5Uniform, dim, gslc, ib)

Convert global slice `gslc` to local slice `lslc` on a given block index `ib`.
"""
function global_slice_to_local_slice(
        file::_BatsrusHDF5Uniform, dim::Int, gslc::OrdinalRange,
        ib::Int
    )
    return trimslice(
        gslc.start - file.common.ncb[dim] * ib, gslc.stop - file.common.ncb[dim] * ib,
        gslc.step, file.common.ncb[dim]
    )
end

"""
Fill `output` by iterating over the block grid defined by `ibx_`, `iby_`, `ibz_`.

`getblock(lslcx, lslcy, lslcz, ib)` returns a block slice for a given block index.
"""
function _fill_output!(
        getblock, output, file::_BatsrusHDF5Uniform,
        gslcx, gslcy, gslcz, ibx_, iby_, ibz_
    )
    nbx, nby, _ = file.nb
    for ibz in ibz_
        lslcz, iz_ = prep_extract_per_block(file, 3, gslcz, ibz)
        for iby in iby_
            lslcy, iy_ = prep_extract_per_block(file, 2, gslcy, iby)
            for ibx in ibx_
                lslcx, ix_ = prep_extract_per_block(file, 1, gslcx, ibx)
                ib = (ibz * nby + iby) * nbx + ibx + 1
                output[ix_, iy_, iz_] = getblock(lslcx, lslcy, lslcz, ib)
            end
        end
    end
    return
end

"""
    extract_var(file::BatsrusHDF5Uniform, var::String; kwargs...)

Extract variable `var` from HDF5 `file`.

# Keywords

  - `xmin`: minimum extracted coordinate in x.
  - `xmax`: maximum extracted coordinate in x.
  - `stepx`: extracted stride in x.
  - `ymin`: minimum extracted coordinate in y.
  - `ymax`: maximum extracted coordinate in y.
  - `stepy`: extracted stride in y.
  - `zmin`: minimum extracted coordinate in z.
  - `zmax`: maximum extracted coordinate in z.
  - `stepz`: extracted stride in z.
  - `verbose::Bool=false`: display type and size information of output variable.
  - `lazy::Bool=false`: if `true`, read data per-block instead of loading the entire dataset.
    Use for files too large to fit in memory.
"""
function Batsrus.extract_var(
        file::_BatsrusHDF5Uniform{TI, TF}, var::String;
        xmin = -Inf32, xmax = Inf32, stepx::Int = 1,
        ymin = -Inf32, ymax = Inf32, stepy::Int = 1,
        zmin = -Inf32, zmax = Inf32, stepz::Int = 1,
        verbose::Bool = false, lazy::Bool = false
    ) where {TI, TF}
    gslcx, xl_new,
        xu_new, ibx_ = prep_extract(file; dim = 1, vmin = xmin, vmax = xmax, step = stepx)
    gslcy, yl_new,
        yu_new, iby_ = prep_extract(file; dim = 2, vmin = ymin, vmax = ymax, step = stepy)
    gslcz, zl_new,
        zu_new, ibz_ = prep_extract(file; dim = 3, vmin = zmin, vmax = zmax, step = stepz)
    nsize = (length(gslcx), length(gslcy), length(gslcz))

    dset = file.common.fid[var]
    output = Array{TF, 3}(undef, nsize)

    if verbose
        @info "output $(typeof(output))"
        @info "nx, ny, nz = $nsize"
    end

    if !lazy
        # Fast path: bulk read + in-memory slicing
        input = read(dset)::Array{TF, 4}
        _fill_output!(output, file, gslcx, gslcy, gslcz, ibx_, iby_, ibz_) do lslcx, lslcy, lslcz, ib
            @view input[lslcx, lslcy, lslcz, ib]
        end
    else
        # Memory-safe path: per-block partial reads
        _fill_output!(output, file, gslcx, gslcy, gslcz, ibx_, iby_, ibz_) do lslcx, lslcy, lslcz, ib
            dset[lslcx, lslcy, lslcz, ib]
        end
    end

    return output, (xl_new, yl_new, zl_new), (xu_new, yu_new, zu_new)
end

end
