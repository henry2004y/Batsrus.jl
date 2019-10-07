# Script for converting BATSRUS unstructured binary outputs to VTK formats.
# Example usage:
# julia -p 8 convert_parallel.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/home1/06426/hyzhou/VisAnaJulia");
@everywhere using VisAna, Glob

filenamesIn = "cut*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
