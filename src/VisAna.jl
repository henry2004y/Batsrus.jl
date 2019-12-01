module VisAna
# Reimplementation of BATSRUS data reader in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, PyPlot, Printf, PyCall, Dierckx, WriteVTK

export readdata, readlogdata, plotdata, plotlogdata, animatedata, readtecdata,
       Data, FileList, convertVTK, get_vars#,
       plot, scatter, contour, contourf, plot_surface, tricontourf,
       plot_trisurf, streamplot

import PyPlot.plot, PyPlot.scatter, PyPlot.contour, PyPlot.contourf,
       PyPlot.plot_surface, PyPlot.tricontourf, PyPlot.plot_trisurf,
       PyPlot.streamplot


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
include("visual.jl")

end
