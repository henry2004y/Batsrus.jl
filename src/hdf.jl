"Module for BATSRUS HDF5 file processing."
module HDF

using HDF5

export BatsrusHDF5Uniform, extract_field, squeeze

abstract type BatsrusHDF5File end

"""
BATSRUS hdf5 file wrapper.

The data are stored in blocks, i.e., each field component is stored in a 4d array in the
order (iblock, iz, iy, ix). This is a generic wrapper and does not assume grid type, i.e.,
uniform, stretched nonuniform, or AMR, etc. Classes to handle data with different grids can
be derived from this class.
"""
struct HDF5Common{TI<:Signed, TF<:AbstractFloat} <: BatsrusHDF5File
   fid::HDF5.File
   version::TI
   "Type of geometry. 0 for Cartesian."
   geometry::TI
   time::TF
   "Minimum coordinates for the whole grid."
   coordmin::Vector{TF}
   "Maximum coordinates for the whole grid."
   coordmax::Vector{TF}
   timestep::TI
   "True dimension of the data despite the stored data dimension."
   ndim::TI
   "Number of cell in each block along x, y, z."
   ncb::Vector{TI}
   "If boundary condition is periodic."
   isperiodic::Vector{Bool}

   multi_cell_dims::Vector{Bool}
   "Lengths along each direction for the whole grid."
   extent::Vector{TF}

   function HDF5Common(filename::AbstractString)
      fid = h5open(filename, "r")
      meta_int = read(fid["Integer Plot Metadata"])::Vector{Int32}
      meta_real = read(fid["Real Plot Metadata"])::Vector{Float32}

      version = meta_int[1]
      geometry = meta_int[11]
      time = meta_real[1]
      coordmin = meta_real[2:2:6]
      coordmax = meta_real[3:2:7]
      timestep = meta_int[2]
      ndim = meta_int[3]
      ncb = meta_int[8:10]
      isperiodic = meta_int[12:14]

      single_cell_dims = ncb .== 0
      multi_cell_dims = .!single_cell_dims
      extent = coordmax .- coordmin

      new{eltype(meta_int), eltype(meta_real)}(fid, version, geometry, time, coordmin,
         coordmax, timestep, ndim, ncb, isperiodic, multi_cell_dims, extent)
   end
end


struct BatsrusHDF5Uniform{TI, TF} <: BatsrusHDF5File
   common::HDF5Common{TI, TF}
   "Numbers of cells along each direction"
   nc::Vector{TI}
   "Number of blocks along each direction"
   nb::Vector{TI}
   "Grid resolution along each direction"
   dcoord::Vector{TF}
   "Block length along each direction"
   dblock::Vector{TF}

   function BatsrusHDF5Uniform(filename::AbstractString)
      bf = HDF5Common(filename)

      nc = Int32[1,1,1]
      try
         extent::Matrix{Int32} = bf.fid["MinLogicalExtents"] |> read
         nc[bf.multi_cell_dims] = extent[:, end] + bf.ncb[bf.multi_cell_dims]
      catch
         iCoord_DB::Matrix{Int32} = bf.fid["iCoord_DB"] |> read
         nc[bf.multi_cell_dims] = iCoord_DB[:,end] + bf.ncb[bf.multi_cell_dims]
      end
      nb = nc .÷ bf.ncb

      dcoord = bf.extent ./ nc
      dblock = bf.extent ./ nb

      TI, TF = findparam(bf)

      new{TI, TF}(bf, nc, nb, dcoord, dblock)
   end
end

findparam(::HDF5Common{TI, TF}) where {TI, TF} = (TI, TF)


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
function prep_extract(file::BatsrusHDF5Uniform;
   dim::Int=1, vmin=-Inf32, vmax=Inf32, step::Int=1)
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

"Return lower corner index."
function coord2index(file::BatsrusHDF5Uniform{TI, TF}, dim::Int, x::Real) where {TI, TF}
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
function prepslice(file::BatsrusHDF5Uniform; dim::Int=1, vmin, vmax, step::Int=1)
   start = coord2index(file, dim, vmin)
   stop = max(coord2index(file, dim, vmax), start)

   slc_new = trimslice(start, stop, step, file.nc[dim])
   xl_new = file.common.coordmin[dim] + (slc_new.start-1) * file.dcoord[dim]
   xu_new = file.common.coordmin[dim] + slc_new.stop * file.dcoord[dim]

   slc_new, xl_new, xu_new
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

   start:step:stop
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
@inline function prep_extract_per_block(file::BatsrusHDF5Uniform, dim::Int,
   gslc::OrdinalRange, ib::Int)
   # compute local slice on this block
   lslc = global_slice_to_local_slice(file, dim, gslc, ib)
   # compute the index bounds in the output array
   ix0 = length(gslc.start:gslc.step:(lslc.start + file.common.ncb[dim]*ib))
   ix1 = ix0 + length(lslc) - 1
   
   lslc, ix0:ix1
end

"""
    global_slice_to_local_slice(file::BatsrusHDF5Uniform, dim, gslc, ib)

Convert global slice `gslc` to local slice `lslc` on a given block index `ib`.
"""
function global_slice_to_local_slice(file::BatsrusHDF5Uniform, dim::Int, gslc::OrdinalRange,
   ib::Int)

   trimslice(gslc.start - file.common.ncb[dim]*ib, gslc.stop - file.common.ncb[dim]*ib,
      gslc.step, file.common.ncb[dim])
end

"""
    extract_field(file::BatsrusHDF5Uniform, var::String; kwargs...)

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
- `verbose::Bool=true`: display type and size information of output field.
"""
function extract_field(file::BatsrusHDF5Uniform, var::String;
   xmin=-Inf32, xmax=Inf32, stepx::Int=1, ymin=-Inf32, ymax=Inf32, stepy::Int=1,
   zmin=-Inf32, zmax=Inf32, stepz::Int=1, verbose::Bool=false)
   nbx, nby, nbz = file.nb

   gslcx, xl_new, xu_new, ibx_ = prep_extract(file; dim=1, vmin=xmin, vmax=xmax, step=stepx)
   gslcy, yl_new, yu_new, iby_ = prep_extract(file; dim=2, vmin=ymin, vmax=ymax, step=stepy)
   gslcz, zl_new, zu_new, ibz_ = prep_extract(file; dim=3, vmin=zmin, vmax=zmax, step=stepz)
   nsize = (length(gslcx), length(gslcy), length(gslcz))

   input = read(file.common.fid[var])::Array{Float32, 4}
   output = Array{eltype(input), 3}(undef, nsize)

   if verbose
      @info "output $(typeof(output))"
      @info "nx, ny, nz = $nsize"
   end

   for ibz in ibz_
      lslcz, iz_ = prep_extract_per_block(file, 3, gslcz, ibz)
      for iby in iby_
         lslcy, iy_ = prep_extract_per_block(file, 2, gslcy, iby)
         for ibx in ibx_
            lslcx, ix_ = prep_extract_per_block(file, 1, gslcx, ibx)
            ib = (ibz * nby + iby) * nbx + ibx + 1
            output[ix_, iy_, iz_] = @view input[lslcx, lslcy, lslcz, ib]
         end
      end
   end

   output, (xl_new, yl_new, zl_new), (xu_new, yu_new, zu_new)
end

"Squeeze singleton dimensions for an array `A`."
function squeeze(A::AbstractArray)
   singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)
   
   dropdims(A, dims=singleton_dims)
end

end