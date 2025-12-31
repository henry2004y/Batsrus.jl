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
   "z=0_fluid_region0_0_t00001640_n00010142.out")

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
SUITE["read"]["Extract density"] = @benchmarkable Batsrus.getvar($bd, "Rho")
SUITE["read"]["Extract Bmag"] = @benchmarkable get_magnitude($bd, :B)
SUITE["read"]["Interp2d"] = @benchmarkable Batsrus.interp2d($bd, "rho")

file = joinpath(directory, files[3])
bd = load(file)
SUITE["read"]["Anisotropy"] = @benchmarkable get_anisotropy($bd, 1)

function generate_mock_amrex_data(output_dir::String)
   ptype = "particles"
   base_dir = joinpath(output_dir, ptype)
   if isdir(base_dir)
      rm(base_dir, recursive = true)
   end
   mkpath(base_dir)

   # Create Header
   header_path = joinpath(base_dir, "Header")
   open(header_path, "w") do f
      println(f, "Version_double")
      println(f, "3") # dim
      println(f, "2") # num_real_extra (total real = 3 + 2 = 5)
      println(f, "u")
      println(f, "v")
      println(f, "2") # num_int_extra (total int = 2 + 2 = 4)
      println(f, "id_1")
      println(f, "id_2")
      println(f, "0") # is_checkpoint (False)
      println(f, "1000") # num_particles (larger for benchmark)
      println(f, "1001") # max_next_id
      println(f, "0") # finest_level
      println(f, "1") # grids_per_level[0]
      println(f, "1 1000 0") # grid info: which, count, where
   end

   # Create Level directory
   level_dir = joinpath(base_dir, "Level_0")
   mkpath(level_dir)

   # Create Particle_H
   particle_h_path = joinpath(level_dir, "Particle_H")
   open(particle_h_path, "w") do f
      println(f, "(1 0") # num_boxes level
      println(f, "((0,0,0) (10,10,10) (0,0,0))")
   end

   # Create Data file
   data_fn = joinpath(level_dir, "DATA_00001")
   open(data_fn, "w") do f
      # Write 1000 particles
      # structure: 1000 * (3+2) reals = 5000 doubles
      data = zeros(Float64, 5000)
      for i in 1:1000
         # x, y, z
         data[(i - 1) * 5 + 1] = Float64(i % 10)
         data[(i - 1) * 5 + 2] = Float64(i % 10)
         data[(i - 1) * 5 + 3] = Float64(i % 10)
         # u, v
         data[(i - 1) * 5 + 4] = Float64(rand())
         data[(i - 1) * 5 + 5] = Float64(rand())
      end
      write(f, data)
   end

   # Create Main Header (for domain info)
   main_header_path = joinpath(output_dir, "Header")
   open(main_header_path, "w") do f
      println(f, "HyperCLaw-V1.1")
      println(f, "0") # num_fields
      println(f, "3") # dim
      println(f, "0.0") # time
      println(f, "0") # refine_ratio
      println(f, "0.0 0.0 0.0") # left_edge
      println(f, "10.0 10.0 10.0") # right_edge
      println(f, "0")
      println(f, "((0,0,0) (10,10,10) (0,0,0))") # domain size
   end
end

println("Generating mock AMReX data for benchmark...")
amrex_dir = joinpath(directory, "amrex_mock")
mkpath(amrex_dir)
generate_mock_amrex_data(amrex_dir)

SUITE["amrex"] = BenchmarkGroup()
SUITE["amrex"]["load"] = @benchmarkable AMReXParticleData($amrex_dir)

# Helper for selection benchmark (pre-load data)
amrex_data = AMReXParticleData(amrex_dir)
Batsrus.load_data!(amrex_data) # Force load for selection benchmark

SUITE["amrex"]["select_region"] = @benchmarkable select_particles_in_region($amrex_data, x_range=(2.5, 4.5))

