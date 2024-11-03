module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using LinearAlgebra: normalize, ×, ⋅
using Printf, Reexport, Requires
using Parsers
using Interpolations: cubic_spline_interpolation, BSpline, Linear, scale, interpolate
import NaturalNeighbours as NN
using StaticArrays: SVector, @SMatrix, SA

export BATLData,
   load, readlogdata, readtecdata, showhead, # io
   getvars, getvar, cutdata, subvolume, subsurface, # select
   Batl, convertTECtoVTU, convertIDLtoVTK, readhead, readtree, getConnectivity, # vtk
   interp1d, interp2d, slice1d, get_var_range, squeeze, get_range # plot/utility

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
struct BATLData{dim, T<:AbstractFloat}
   "header information"
   head::NamedTuple
   "grid"
   x::Array{T}
   "variables"
   w::Array{T}
   "file information"
   list::FileList
end

include("unit/UnitfulBatsrus.jl")
using .UnitfulBatsrus

include("io.jl")
include("select.jl")
include("vtk.jl")
include("utility.jl")
include("plot/plots.jl")

include("hdf.jl")
@reexport using .HDF

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
end

end
