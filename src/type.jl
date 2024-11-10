abstract type AbstractBATS{Dim, T} end

@enum FileType Real4Bat=1 Real8Bat=2 AsciiBat=3 LogBat=4 TecBat=5

"Type for Batsrus file information."
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

"Batsrus file head information."
struct BatsHead
   ndim::Int32
   headline::SubString{String}
   it::Int32
   time::Float32
   gencoord::Bool
   neqpar::Int32
   nw::Int32
   nx::Vector{Int32}
   eqpar::Vector{Float32}
   variables::Vector{SubString{String}}
   wnames::Vector{SubString{String}}
end

"Batsrus data container."
struct BATS{Dim, Dimp1, T<:AbstractFloat} <: AbstractBATS{Dim, T}
   head::BatsHead
   x::Array{T, Dimp1}
   w::Array{T, Dimp1}
   list::FileList

   function BATS(head, x::Array{T, Dimp1}, w::Array{T, Dimp1}, list) where {T, Dimp1}
      @assert Int(head.ndim) + 1 == Dimp1 "Dimension mismatch!"
      new{Int(head.ndim), Dimp1, T}(head, x, w, list)
   end
end