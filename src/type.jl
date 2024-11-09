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

@computed struct BATS{Dim, T<:AbstractFloat} <: AbstractBATS{Dim, T}
   head::BatsHead
   x::Array{T, Dim+1}
   w::Array{T, Dim+1}
   list::FileList
end