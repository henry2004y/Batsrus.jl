module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Printf, Requires

export Data

"""
	FileList

Type for the file information. Contains filename `name`, file type `type` of
ascii or binary, file size `bytes`, and number of snapshots `npictinfiles`.
"""
struct FileList
   name::String
   type::String
   bytes::Int64
   npictinfiles::Int64
end

"""
    Data

Primary data storage type, with fields `head` of header info, grid `x`, value
`w`, and file info `list`.
"""
struct Data{T<:AbstractFloat}
   head::NamedTuple
   x::Array{T}
   w::Array{T}
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
