module VisAna
# Reimplementation of BATSRUS data reader in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Printf, PyCall, Dierckx, WriteVTK

export readdata, readlogdata, readtecdata, Data, FileList, convertVTK,
       plot

struct Data{T}
   x::Array{T}
   w::Array{T}
end

struct FileList
   name::String
   type::String
   bytes::Int64
   npictinfiles::Int64
end

struct Vars
   data::Dict{String, Array{Float32}}
end

include("io.jl")
#include("visual.jl")
include("visual_plot.jl")

end
