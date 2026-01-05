# Tests of BATSRUS.jl

using Batsrus, Test, SHA, LazyArtifacts
using Batsrus.UnitfulBatsrus, Unitful
using Batsrus: At, Near # DimensionalData
using RecipesBase
using Suppressor: @capture_out, @capture_err, @suppress_out, @suppress_err, @suppress
using CairoMakie
using CairoMakie
# Check if we should run PyPlot tests (Only on Linux in CI, or always locally)
const RUN_PYPLOT_TESTS = Sys.islinux() || get(ENV, "CI", "false") != "true"
if RUN_PYPLOT_TESTS
   using PyPlot
end
using FHist

ENV["MPLBACKEND"] = "agg" # no GUI

function filecmp(path1::AbstractString, path2::AbstractString)
   stat1, stat2 = stat(path1), stat(path2)
   if !(isfile(stat1) && isfile(stat2)) || filesize(stat1) != filesize(stat2)
      return false
   end
   stat1 == stat2 && return true # same file
   open(path1, "r") do file1
      open(path2, "r") do file2
         buf1 = Vector{UInt8}(undef, 32768)
         buf2 = similar(buf1)
         while !eof(file1) && !eof(file2)
            n1 = readbytes!(file1, buf1)
            n2 = readbytes!(file2, buf2)
            n1 != n2 && return false
            0 != Base._memcmp(buf1, buf2, n1) && return false
         end
         return eof(file1) == eof(file2)
      end
   end
end

datapath = artifact"testdata"

@testset "Batsrus.jl" begin
   tests = isempty(ARGS) ? ["basic", "io", "plotting", "analysis"] : ARGS
   for t in tests
      file = joinpath(@__DIR__, "tests_$(t).jl")
      if isfile(file)
         include(file)
      else
         @warn "Test file $file not found!"
      end
   end
end
