module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Printf, Requires

export Data

"Type for the file information."
struct FileList
   "filename"
   name::String
   "file type"
   type::String
   "file size"
   bytes::Int
   "number of snapshots"
   npictinfiles::Int
end

"Primary data storage type"
struct Data{T<:AbstractFloat}
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

include("unit/UnitfulBatsrus.jl")
using UnitfulRecipes, .UnitfulBatsrus

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
   @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
      include("plot/plots.jl")
   end
end

end
