# Tests of BATSRUS.jl

using Batsrus, Test, SHA, LazyArtifacts
using Batsrus.UnitfulBatsrus, Unitful
using Batsrus: At, Near # DimensionalData
using RecipesBase
using Suppressor: @capture_out, @capture_err, @suppress_out, @suppress_err
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

@testset "Batsrus.jl" begin
   datapath = artifact"testdata"
   @testset "Reading 1D ascii" begin
      file = "1d__raw_2_t25.60000_n00000258.out"
      bd = @suppress_err begin
         load(joinpath(datapath, file), verbose = true)
      end
      @test startswith(repr(bd), "filename : 1d")
      @test Batsrus.setunits(bd.head, "NORMALIZED")
      @test extrema(bd.x) == (-127.5, 127.5)
      @test extrema(bd.w) == (-0.79960780498, 1.9394335293)
   end

   @testset "Reading 2D structured binary" begin
      file = "z=0_raw_1_t25.60000_n00000258.out"
      @test_throws ArgumentError load(joinpath(datapath, file), npict = 2)
      bd = load(joinpath(datapath, file))
      @test bd.head.time == 25.6f0
      @test extrema(bd.x) == (-127.5f0, 127.5f0)
      @test extrema(bd.w) == (-0.79985905f0, 1.9399388f0)
      plotrange = [-10.0, 10.0, -Inf, Inf]
      x, y, w = interp2d(bd, "rho", plotrange)
      @test w[1, end] == 0.6848635077476501
      @test Batsrus.fill_vector_from_scalars(bd, :B)[:, end, end] ==
            Float32[1.118034, -0.559017, 0.0]
      @test get_magnitude(bd, :B)[128, 2] == 0.9223745f0
      @test get_magnitude2(bd, :B)[128, 2] == 0.8507747f0
      # Linear interpolation at a given point
      d = interp1d(bd, "rho", Float32[0.0, 0.0])
      @test d == 0.6936918f0
      # Linear interpolation along a line
      point1 = Float32[-10.0, -1.0]
      point2 = Float32[10.0, 1.0]
      w = interp1d(bd, "rho", point1, point2)
      @test sum(w) == 14.942542f0
      w = slice1d(bd, "rho", 1, 1)
      @test sum(w) == 4.0f0
      @test get_var_range(bd, "rho") == (0.11626893f0, 1.0f0)
      @test length(bd["rho"][X = -10 .. 10, Y = Near(0.9)]) == 20
      @test bd["rho"][X = At(0.5), Y = At(0.5)] == 0.6940014f0

      file = "z=0_fluid_region0_0_t00001640_n00010142.out"
      bd = load(joinpath(datapath, file))
      x, y = Batsrus.meshgrid(bd)
      @test length(x) == 601 && y[2] == 0.0f0
      x, y = Batsrus.meshgrid(bd, Float32[-100, 100, -Inf, Inf])
      @test length(x) == 4
      @test get_magnitude(bd, :E)[2, 1] == 2655.4805f0
      @test get_magnitude2(bd, :E)[2, 1] == 7.051577f6
      @test Batsrus.fill_vector_from_scalars(bd, :E)[:, 2, 1] ==
            Float32[-241.05942, -2644.2058, -40.53219]
      @test get_magnitude2(bd, :U0)[2, 1] == 33784.973f0
      anisotropy_s0 = get_anisotropy(bd, 0)[1:2, 1]
      @test anisotropy_s0 ≈ Float32[1.2630985, 2.4700143]
      @test get_anisotropy(bd, 0, method = :rotation)[1:2, 1] ≈ anisotropy_s0
      @test get_anisotropy(bd, 1)[1:2, 1] ≈ Float32[1.2906302, 2.6070855]
      w = get_convection_E(bd)
      @test w[2][2, 1] ≈ -2454.3933f0
      w = get_hall_E(bd)
      @test w[2][2, 1] ≈ -782.2945f0
      w = get_timeseries([joinpath(datapath, file)], [0.0, 0.0])
      @test w[2][end] == 17.973747f0
   end

   @testset "Reading 2D unstructured ascii" begin
      file = "bx0_mhd_6_t00000100_n00000352.out"
      bd = load(joinpath(datapath, file))
      plotrange = [-Inf, Inf, -Inf, Inf]
      x, y = Batsrus.meshgrid(bd)
      @test length(x) == 117 && length(y) == 246
      x, y, w = interp2d(bd, "rho", plotrange, useMatplotlib = false)
      @test w[1, 2] == 5.000018304080387
      @test get_magnitude(bd, :U)[2] == 71.85452748407637
      @test get_magnitude2(bd, :U)[2] == 5163.073119959886
   end

   @testset "Reading 3D structured binary" begin
      file = "3d_raw.out"
      bd = load(joinpath(datapath, file))
      plotrange = [-50.0, 50.0, -0.5, 0.5]
      X, Z, p = cutdata(bd, "p"; dir = "y", sequence = 1, plotrange)
      @test p[1] ≈ 0.560976f0 && p[2] ≈ 0.53704995f0
      @test size(bd["p"]) == (8, 8, 8)
   end
   #TODO: add tecplot tests
   #@testset "Reading Tecplot" begin

   #end

   @testset "Log" begin
      logfilename = joinpath(datapath, "log_n000001.log")
      head, bd = readlogdata(logfilename)
      @test isa(head, NamedTuple)
      @test extrema(bd) == (-0.105, 258.0)
   end

   @testset "VTK" begin
      file = joinpath(datapath, "z=0_fluid_region0_0_t00001640_n00010142.out")
      convertIDLtoVTK(file)
      sha_str = bytes2hex(open(sha1, "out.vti"))
      @test sha_str == "7c7540e684aaf1823b7e2792f555e89ad98bcd31"

      file = joinpath(datapath, "3d_raw.out")
      convertIDLtoVTK(file)
      sha_str = bytes2hex(open(sha1, "out.vti"))
      @test sha_str == "8d535c0656daa9224aab048d036da0cc2f667bdc"
      rm("out.vti")

      file = joinpath(datapath, "3d_bin.dat")
      convertTECtoVTU(file)
      sha_str = bytes2hex(open(sha1, "out.vtu"))
      @test sha_str == "5b04747666542d802357dec183177f757754a254"
      rm("out.vtu")

      filetag = joinpath(datapath, "3d_mhd_amr/3d__mhd_1_t00000000_n00000000")
      batl = Batl(readhead(filetag * ".info"), readtree(filetag)...)
      # local block index check
      @test Batsrus.find_grid_block(batl, [1.0, 0.0, 0.0]) == 4 &&
            Batsrus.find_grid_block(batl, [100.0, 0.0, 0.0]) == -100
      connectivity = getConnectivity(batl)
      sha_str = bytes2hex(sha256(string(connectivity)))
      @test sha_str == "c6c5a65a46d86a9ba4096228c1516f89275e45e295cd305eb70c281a770ede74"
   end

   @testset "HDF5" begin
      file = "3d__var_3_n00000000_single.batl"
      bd = BatsrusHDF5Uniform(joinpath(datapath, file))
      var, _, _ = extract_var(bd, "bx")
      @test size(var) == (8, 8, 4) && eltype(var) == Float32
      file = "3d__var_3_n00000000_double.batl"
      bd = BatsrusHDF5Uniform(joinpath(datapath, file))
      var, _, _ = extract_var(bd, "bx")
      @test size(var) == (8, 8, 4) && eltype(var) == Float64
   end

   @testset "Plotting" begin
      @testset "Plots" begin
         RecipesBase.is_key_supported(k::Symbol) = true
         # 1D
         file = "1d__raw_2_t25.60000_n00000258.out"
         bd = load(joinpath(datapath, file))
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bd, "Rho")
         @test getfield(rec[1], 1)[:seriestype] == :path

         file = "z=0_raw_1_t25.60000_n00000258.out"
         bd = load(joinpath(datapath, file))
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bd, "p")
         @test getfield(rec[1], 1)[:seriestype] == :contourf
      end

      @testset "Makie" begin
         file = "1d__raw_2_t25.60000_n00000258.out"
         bd = load(joinpath(datapath, file))
         fig, ax, plt = CairoMakie.lines(bd, "Rho")
         @test plt isa Lines

         file = "z=0_raw_1_t25.60000_n00000258.out"
         bd = load(joinpath(datapath, file))
         fig, ax, plt = CairoMakie.heatmap(bd, "p")
         @test plt isa Heatmap
      end

      if RUN_PYPLOT_TESTS
         @testset "PyPlot" begin
            @test size(squeeze(zeros(2, 3, 1))) == (2, 3)
            # 1D ascii
            file = "1d__raw_2_t25.60000_n00000258.out"
            bd = load(joinpath(datapath, file), verbose = false)
            c = PyPlot.plot(bd, "p")
            @test c[1].get_xdata() ≈ bd.x
            @test c[1].get_ydata() ≈ bd.w[:, 10]

            # 2D structured binary
            file = "z=0_raw_1_t25.60000_n00000258.out"
            bd = load(joinpath(datapath, file))
            c = PyPlot.streamplot(bd, "bx;by")
            @test c.lines.get_segments()[2][3] ≈ -118.68871477694084
            c = PyPlot.contourf(bd, "p")
            @test c.get_array()[end] == 0.9750000000000002
            c = @suppress_err PyPlot.contourf(bd, "rho", innermask = true)
            @test c.get_array()[end] == 0.9750000000000002
            c = PyPlot.contour(bd, "rho")
            @test c.get_array()[end] == 1.0500000000000003
            c = PyPlot.contour(bd, "rho"; levels = [1.0])
            @test c.get_array()[end] == 1.0
            c = PyPlot.tricontourf(bd, "rho")
            @test c.get_array()[end] == 0.9750000000000002
            PyPlot.tripcolor(bd, "rho")
            @test isa(gca(), PyPlot.PyObject)
            p = PyPlot.pcolormesh(bd, "p").get_array()
            @test p[end] == 0.1f0
            p = PyPlot.imshow(bd, "p").get_array()
            @test p[2, 128] == 0.51229393f0
            plt.close()
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
            plot_surface(bd, "rho")
            @test isa(gca(), PyPlot.PyObject)
            plt.close()

            # 2D AMR Cartesian
            file = "bx0_mhd_6_t00000100_n00000352.out"
            bd = load(joinpath(datapath, file))
            pcolormesh(bd, "P")
            @test isa(gca(), PyPlot.PyObject)
         end
      end
   end

   @testset "Units" begin
      @test 1.0bu"R" > 2.0bu"Rg"
      file = "y=0_var_1_t00000000_n00000000.out"
      bd = load(joinpath(datapath, file))
      varunit = getunit(bd, "Rho")
      @test varunit == bu"amucc"
      varunit = getunit(bd, "Ux")
      @test varunit == u"km/s"
      varunit = getunit(bd, "Bx")
      @test varunit == u"nT"
      varunit = getunit(bd, "P")
      @test varunit == u"nPa"
      varunit = getunit(bd, "jx")
      @test dimension(varunit) == dimension(Unitful.A / Unitful.m^2)
      varunit = getunit(bd, "ex")
      @test varunit == u"mV/m"
   end

   @testset "AMReX Loader" begin
      mktempdir() do tmpdir
         Batsrus.generate_mock_amrex_data(tmpdir)

         # Test Loading
         data = AMReXParticle(tmpdir)
         @test data.dim == 3
         @test data.header.num_levels == 1
         @test data.header.num_particles == 10

         # Test Data Access
         rdata = data.rdata
         @test size(rdata) == (5, 10)
         @test rdata[1, 1] == 1.0 # x of first particle
         @test rdata[4, 1] == 10.0 # u of first particle
         @test rdata[5, 10] == 1000.0 # v of last particle

         # Test Region Selection
         selected = select_particles_in_region(data, x_range = (2.5, 4.5))
         @test size(selected, 2) == 2
         @test selected[1, 1] == 3.0
         @test selected[1, 2] == 4.0

         # Test Plotting Helper
         h = get_phase_space_density(data, "x", "u", bins = 2)
         @test size(h.bincounts) == (2, 2)

         # Test Weighted Plotting
         Batsrus.generate_mock_amrex_data(tmpdir;
            num_particles = 10,
            real_component_names = ["u", "v", "weight"],
            particle_gen = (i, n_reals) -> (
               Float64(i), 10.0, 10.0, # x, y, z
               Float64(i), 0.0, # u, v
               Float64(2.0) # weight
            )
         )
         data_w = AMReXParticle(tmpdir)
         h_w = get_phase_space_density(data_w, "x", "u", bins = 1,
            x_range = (0.0, 11.0), y_range = (0.0, 11.0))

         # 10 particles, each weight 2.0. Total sum should be 20.0.
         # They all fall into the same bin since we only have 1 bin covering the range.
         @test sum(h_w.bincounts) ≈ 20.0
      end
   end

   @testset "Particle Classification" begin
      mktempdir() do tmpdir
         n_core = 1000
         n_halo = 100
         n_total = n_core + n_halo

         Batsrus.generate_mock_amrex_data(tmpdir;
            num_particles = n_total,
            real_component_names = ["vx", "vy", "vz"],
            particle_gen = (i, n_reals) -> begin
               # Positions (random in box 0-2 (since mock data is small?))
               pos = rand(3) .* 2.0 .- 0.5

               # Velocities
               if i <= n_core
                  # Core velocities (Gaussian)
                  vel = randn(3)
               else
                  # Halo velocities (Broad Gaussian)
                  vel = randn(3) .* 5.0
               end

               (pos..., vel...)
            end
         )

         data = AMReXParticle(tmpdir)
         # Ensure data is loaded
         Batsrus.load_data!(data)

         # Test 1: Classification with known params
         # vth=1.0, u=0
         # nsigma=3.0

         core, halo = classify_particles(
            data; vdim = 3, bulk_vel = [0.0, 0.0, 0.0], vth = 1.0, nsigma = 3.0)

         # Core efficiency: erte of 3 sigma is 0.997.
         # So expected core count ~ n_core * 0.997 + halo_leak
         # Halo leak: halo sigma=5. threshold=3. fraction within 3/5=0.6 sigma of halo distribution?
         # erf(0.6/sqrt(2)) is small.

         # Just check counts roughly
         @test size(core, 2) > n_core * 0.9
         @test size(halo, 2) < n_halo + n_core * 0.1
         @test size(halo, 2) > n_halo * 0.5 # Halo shouldn't be empty

         # Test 2: Auto validation of bulk velocity
         core_auto, halo_auto = classify_particles(
            data; vdim = 3, bulk_vel = nothing, vth = 1.0, nsigma = 3.0)

         @test size(core_auto, 2)≈size(core, 2) atol=50

         # Test 3: 1D
         core_1d, halo_1d = classify_particles(
            data; vdim = 1, bulk_vel = [0.0], vth = 1.0, nsigma = 3.0)

         # In 1D, we verify against x-velocity only.
         # Similar statistics hold.
         @test size(core_1d, 2) > n_core * 0.9
      end
   end

   @testset "Field Aligned Transform" begin
      # Test with simple magnetic field along X
      # v_para should be vx, v_perp should be sqrt(vy^2 + vz^2)
      b_field = [1.0, 0.0, 0.0]
      transform = Batsrus.get_particle_field_aligned_transform(b_field)

      # Mock data: 2 particles, 6 components (x, y, z, vx, vy, vz)
      # Particle 1: vx=1, vy=0, vz=0 => v_para=1, v_perp=0
      # Particle 2: vx=0, vy=3, vz=4 => v_para=0, v_perp=5
      names = ["x", "y", "z", "vx", "vy", "vz"]
      data = [0.0 0.0 0.0 1.0 0.0 0.0;
              0.0 0.0 0.0 0.0 3.0 4.0]'

      new_data, new_names = transform(data, names)

      @test new_names == ["v_parallel", "v_perp"]
      @test new_data[1, 1] ≈ 1.0
      @test new_data[1, 2] ≈ 0.0
      @test new_data[2, 1] ≈ 0.0
      @test new_data[2, 2] ≈ 5.0

      # Test with E and B
      # B = x, E = y
      # ExB = z
      # b_hat = (1,0,0)
      # d_hat = (0,0,1)
      # e_hat = d x b = (0,0,1) x (1,0,0) = (0,1,0) (which is E direction)
      # So basis is (x, y, z) corresponding to (v_B, v_E, v_BxE)

      b_field = [10.0, 0.0, 0.0]
      e_field = [0.0, 5.0, 0.0]
      transform_eb = Batsrus.get_particle_field_aligned_transform(b_field, e_field)

      # Particle: vx=1, vy=2, vz=3
      # v_B = 1
      # v_E = 2
      # v_BxE = 3
      data_eb = [0.0 0.0 0.0 1.0 2.0 3.0]'

      new_data_eb, new_names_eb = transform_eb(data_eb, names)

      @test new_names_eb == ["v_B", "v_E", "v_BxE"]
      @test new_data_eb[1, 1] ≈ 1.0
      @test new_data_eb[2, 1] ≈ 2.0
      @test new_data_eb[3, 1] ≈ 3.0
   end
end
