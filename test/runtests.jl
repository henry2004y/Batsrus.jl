# Tests of BATSRUS.jl

using Batsrus, Test, SHA, LazyArtifacts
using Batsrus.UnitfulBatsrus, Unitful
using Batsrus: At, Near # DimensionalData
using RecipesBase
using Suppressor: @capture_out, @capture_err, @suppress_out, @suppress_err, @suppress
# Check if we should run PyPlot tests (Linux/Windows CI, or default locally)
const RUN_PYPLOT_TESTS = get(ENV, "TEST_PYPLOT", "true") == "true" &&
                         (Sys.islinux() || Sys.iswindows() ||
                          get(ENV, "CI", "false") != "true")

# Check if we should run Makie tests (MacOS CI, or explicit request)
const RUN_MAKIE_TESTS = get(ENV, "TEST_MAKIE", "false") == "true" || Sys.isapple()

if RUN_MAKIE_TESTS
   using CairoMakie
end

if RUN_PYPLOT_TESTS
   using DimensionalData # for .. (IntervalSet) usage
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
      @testset "$t" begin
         file = joinpath(@__DIR__, "tests_$(t).jl")
         if isfile(file)
            include(file)
         else
            @warn "Test file $file not found!"
         end
      end
   end
end
