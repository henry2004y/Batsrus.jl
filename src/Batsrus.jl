module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using LinearAlgebra: normalize, ×, ⋅, Adjoint
using Printf, Reexport, Requires
using Parsers
using Interpolations: cubic_spline_interpolation, BSpline, Linear, scale, interpolate
import NaturalNeighbours as NN
using StaticArrays: SVector, @SMatrix, SA, MVector
using DimensionalData

export BATS,
   load, readlogdata, readtecdata, showhead, # io
   getvar, cutdata, subvolume, subsurface, get_convection_E, get_hall_E, 
   get_anisotropy, get_vectors, get_magnitude, get_magnitude2,
   fill_vector_from_scalars, # select
   Batl, convertTECtoVTU, convertIDLtoVTK, readhead, readtree, getConnectivity, # vtk
   interp1d, interp2d, slice1d, get_var_range, squeeze, get_range # plot/utility

include("type.jl")
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
