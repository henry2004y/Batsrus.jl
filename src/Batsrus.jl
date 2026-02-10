module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using LinearAlgebra: normalize, ×, ⋅, Adjoint, norm, diag, diagm, tr, dot
using Printf, Reexport
using Parsers
using Interpolations: cubic_spline_interpolation, BSpline, Linear, scale, interpolate
import NaturalNeighbours as NN
using StaticArrays: SVector, @SMatrix, SA, MVector
using DimensionalData
using ProgressMeter

export BATS, BatsrusIDL, BatsrusIDLStructured, BatsrusIDLUnstructured,
    load, readlogdata, readtecdata, showhead, # io
    getvar, cutdata, subvolume, subsurface, get_convection_E, get_hall_E, get_timeseries,
    get_anisotropy, get_vectors, get_magnitude, get_magnitude2,
    fill_vector_from_scalars, # select
    Batl, convertTECtoVTU, convertIDLtoVTK, create_pvd, readhead, readtree, getConnectivity, # vtk
    interp1d, interp2d, slice1d, get_range, get_var_range, squeeze,
    generate_mock_amrex_data, # plot/utility,
    AMReXParticle, AMReXParticleHeader, read_amrex_binary_particle_file,
    select_particles_in_region, get_phase_space_density, classify_particles,
    fit_particle_velocity_gmm, get_particle_field_aligned_transform,
    get_gmm_thermal_velocity, # amrex
    BatsrusHDF5File, BatsrusHDF5Uniform, extract_var # hdf5

include("type.jl")
include("unit/UnitfulBatsrus.jl")
using .UnitfulBatsrus

include("io.jl")
include("select.jl")
include("vtk.jl")
include("utility.jl")
include("amrex.jl")
include("plot/plots.jl")

function BatsrusHDF5Uniform end
function extract_var end

include("precompile.jl")

include("plot/stubs.jl")

end
