# Tests of BATSRUS.jl

using Batsrus, Test, SHA, LazyArtifacts
using Batsrus.UnitfulBatsrus, Unitful
using RecipesBase
using Suppressor: @capture_out, @capture_err, @suppress_out, @suppress_err
using CairoMakie
using PyPlot
ENV["MPLBACKEND"]="agg" # no GUI

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
      @test get_anisotropy(bd, 0)[1:2, 1] ≈ Float32[1.2630985, 2.4700143]
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
      batl = Batl(readhead(filetag*".info"), readtree(filetag)...)
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
      @test dimension(varunit) == dimension(Unitful.A/Unitful.m^2)
      varunit = getunit(bd, "ex")
      @test varunit == u"mV/m"
   end
end
