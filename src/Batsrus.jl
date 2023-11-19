module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Printf, Requires

export BATLData

"Type for the file information."
struct FileList
   "filename"
   name::String
   "file type"
   type::Symbol
   "directory"
   dir::String
   "file size"
   bytes::Int
   "number of snapshots"
   npictinfiles::Int
   "length of meta data"
   lenhead::Int
end

"Primary Batsrus data storage type."
struct BATLData{T<:AbstractFloat}
   "header information"
   head::NamedTuple
   "grid"
   x::Array{T}
   "variables"
   w::Array{T}
   "file information"
   list::FileList
end

include("io.jl")
include("select.jl")
include("vtk.jl")
include("plot/utility.jl")
include("plot/plots.jl")

include("unit/UnitfulBatsrus.jl")
using .UnitfulBatsrus

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
end

end
