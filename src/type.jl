abstract type AbstractBATS{Dim, T} end

@enum FileType Real4Bat=1 Real8Bat=2 AsciiBat=3 LogBat=4 TecBat=5

"""
Type for Batsrus file information.
"""
struct FileList
   "filename"
   name::String
   "file type"
   type::FileType
   "directory"
   dir::String
   "file size"
   bytes::Int
   "number of snapshots"
   npictinfiles::Int
   "length of meta data"
   lenhead::Int
end

"""
Batsrus file head information.
"""
struct BatsHead
   ndim::Int
   headline::SubString{String}
   it::Int
   time::Float32
   gencoord::Bool
   neqpar::Int
   nw::Int
   nx::Vector{Int}
   eqpar::Vector{Float32}
   coord::Vector{SubString{String}}
   wname::Vector{SubString{String}}
   param::Vector{SubString{String}}
end

"Cartesian grid data."
struct BATSCartesian{Dim, TV <: AbstractFloat, TX, TW} <: AbstractBATS{Dim, TV}
   head::BatsHead
   list::FileList
   "Grid"
   x::TX
   "Variables"
   w::TW

   function BATSCartesian(head, list, x::Array{TV, Dimp1}, w::Array{TV, Dimp1}) where {TV, Dimp1}
      @assert head.ndim + 1 == Dimp1 "Dimension mismatch!"
      if head.ndim == 2
         xrange, yrange = get_range(x)
         x = DimArray(x, (X(xrange), Y(yrange), :dim))
         w = DimArray(w, (X(xrange), Y(yrange), Dim{:var}(head.wname)))
      elseif head.ndim == 3
         xrange, yrange, zrange = get_range(x)
         x = DimArray(x, (X(xrange), Y(yrange), Z(zrange), :dim))
         w = DimArray(w, (X(xrange), Y(yrange), Z(zrange), Dim{:var}(head.wname)))
      elseif head.ndim == 1
         (xrange,) = get_range(x)
         x = DimArray(x, (X(xrange), :dim))
         w = DimArray(w, (X(xrange), Dim{:var}(head.wname)))
      end

      new{head.ndim, TV, typeof(x), typeof(w)}(head, list, x, w)
   end
end

"Unstructured grid data."
struct BATSUnstructured{Dim, TV <: AbstractFloat, TX, TW} <: AbstractBATS{Dim, TV}
   head::BatsHead
   list::FileList
   "Grid"
   x::TX
   "Variables"
   w::TW

   function BATSUnstructured(head, list, x::Array{TV, Dimp1}, w::Array{TV, Dimp1}) where {TV, Dimp1}
      @assert head.ndim + 1 == Dimp1 "Dimension mismatch!"
      if head.ndim == 2
         x = DimArray(x, (X, Y, :dim))
         w = DimArray(w, (X, Y, Dim{:var}(head.wname)))
      elseif head.ndim == 3
         x = DimArray(x, (X, Y, Z, :dim))
         w = DimArray(w, (X, Y, Z, Dim{:var}(head.wname)))
      elseif head.ndim == 1
         x = DimArray(x, (X, :dim))
         w = DimArray(w, (X, Dim{:var}(head.wname)))
      end

      new{head.ndim, TV, typeof(x), typeof(w)}(head, list, x, w)
   end
end
