using BenchmarkTools, Downloads, Tar, CodecZlib

t = @elapsed using Batsrus
println("Julia version is $VERSION")
println(string("Batsrus.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Batsrus.jl...")
println()

testdata_url = "https://github.com/henry2004y/batsrus_data/raw/master/batsrus_data.tar.gz"
directory = "data"
files = ("1d__raw_2_t25.60000_n00000258.out", "z=0_raw_1_t25.60000_n00000258.out",
   "z=0_fluid_region0_0_t00001640_n00010142.out", "3d_raw.out")

# Check if all files already exist
if joinpath(directory, files[1]) |> isfile
   println("✅ All data files already exist.")
else
   println("⬇️ Downloading and extracting data...")

   # Download and extract the data
   testdata = Downloads.download(testdata_url)
   open(GzipDecompressorStream, testdata) do io
      Tar.extract(io, directory)
   end

   println("📦 Extraction complete.")
end

const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["IO"])
file = joinpath(directory, files[1])
SUITE["read"]["ASCII"] = @benchmarkable load($file)

file = joinpath(directory, files[2])
bd = load(file)
SUITE["read"]["Load binary structured"] = @benchmarkable load($file)
SUITE["read"]["Extract density"] = @benchmarkable Batsrus.getvar($bd, "rho")
SUITE["read"]["Extract Bmag"] = @benchmarkable get_magnitude($bd, :B)
SUITE["read"]["Interp2d"] = @benchmarkable Batsrus.interp2d($bd, "rho")

file = joinpath(directory, files[3])
bd = load(file)
SUITE["read"]["Anisotropy"] = @benchmarkable get_anisotropy($bd, 1)

file = joinpath(directory, files[4])
bd = load(file)
SUITE["read"]["Cutdir"] = @benchmarkable cutdata($bd, "p", dir = "y", sequence = 1)
SUITE["read"]["Cutdir subset"] = @benchmarkable cutdata(
   $bd, "p", dir = "y", sequence = 1, plotrange = [-50.0, 50.0, -0.5, 0.5])

println("Generating mock AMReX data for benchmark...")
amrex_dir = joinpath(directory, "amrex_mock")
mkpath(amrex_dir)
Batsrus.generate_mock_amrex_data(
   amrex_dir,
   num_particles = 1000,
   real_component_names = ["u", "v", "w", "weight"],
   particle_gen = (i, n_reals) -> (
      Float64(i % 10), Float64(i % 10), Float64(i % 10), # x, y, z
      randn(), randn(), randn(), # u, v, w
      rand() # weight
   )
)

SUITE["amrex"] = BenchmarkGroup()
SUITE["amrex"]["load"] = @benchmarkable AMReXParticle($amrex_dir)

# Helper for selection benchmark (pre-load data)
amrex_data = AMReXParticle(amrex_dir)
Batsrus.load_data!(amrex_data) # Force load for selection benchmark

SUITE["amrex"]["select_region"] = @benchmarkable select_particles_in_region(
   $amrex_data, x_range = (2.5, 4.5))

SUITE["amrex"]["phase_space_3d"] = @benchmarkable get_phase_space_density(
   $amrex_data, "u", "v", "w")

# Helper for file-based selection benchmark (lazy load)
amrex_data_unloaded = AMReXParticle(amrex_dir)
SUITE["amrex"]["select_region_from_files"] = @benchmarkable select_particles_in_region(
   $amrex_data_unloaded, x_range = (2.5, 4.5))
