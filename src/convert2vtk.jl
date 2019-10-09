# Script for converting BATSRUS unstructured binary outputs to VTK formats.
# Example usage:
# julia -p 8 convert_parallel.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/home1/06426/hyzhou/VisAnaJulia");
@everywhere using VisAna, Glob

catcommand = `./pTEC`
run(catcommand)

filenamesIn = "cut*.dat"
dir = "GM"
originDir = pwd()
cd(dir)
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn)
filenames = vcat(filenames, filesfound)
# Do not work on files that have already been converted
filenames = [fname for fname in filenames if ~isfile(fname[1:end-3]*"vtu")]

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end

cd(originDir) # Return to the starting directory

# Choose whether or not to delete the original *.dat files.
