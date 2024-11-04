using BenchmarkTools, LazyArtifacts

t = @elapsed using Batsrus
println("Julia version is $VERSION")
println(string("Batsrus.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Batsrus.jl...")
println()

directory = artifact"testdata"
files = ("1d__raw_2_t25.60000_n00000258.out", "z=0_raw_1_t25.60000_n00000258.out",
   "bx0_mhd_6_t00000100_n00000352.out")


const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["IO"])
file = joinpath(directory, files[1])
SUITE["read"]["ASCII"] = @benchmarkable load($file)

file = joinpath(directory, files[2])
bd = load(file)
SUITE["read"]["Extract density"] = @benchmarkable Batsrus.getvar($bd, "Rho")
SUITE["read"]["Extract B"] = @benchmarkable Batsrus.getvar($bd, "B")
SUITE["read"]["Binary structured"] = @benchmarkable load($file)

file = joinpath(directory, files[3])
SUITE["read"]["Binary unstructured"] = @benchmarkable load($file)
