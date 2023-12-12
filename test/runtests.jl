# Tests of BATSRUS.jl

using Batsrus, Test, SHA, LazyArtifacts
using Batsrus.UnitfulBatsrus, Unitful
using RecipesBase
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
      bd = load(file, dir=datapath, verbose=true)
      @test startswith(repr(bd), "filename : 1d")
      @test Batsrus.setunits(bd.head, "NORMALIZED")
      @test isa(bd.head, NamedTuple)
      @test extrema(bd.x) == (-127.5, 127.5)
      @test extrema(bd.w) == (-0.79960780498, 1.9394335293)
   end

   @testset "Reading 2D structured binary" begin
      file = "z=0_raw_1_t25.60000_n00000258.out"
      bd = load(file, dir=datapath)
      @test bd.head.time == 25.6f0
      @test extrema(bd.x) == (-127.5f0, 127.5f0)
      @test extrema(bd.w) == (-0.79985905f0, 1.9399388f0)
      plotrange = [-10.0, 10.0, -Inf, Inf]
      x, y, w = Batsrus.getdata2d(bd, "rho", plotrange)
      @test w[1,end] == 0.6848978549242021
   end

   @testset "Reading 2D unstructured binary" begin
      #file = "z=0_raw_1_t25.60000_n00000258.out"
      #bd = load(file)
      #TODO test getdata2d on gencoord
   end

   @testset "Reading 3D structured binary" begin
      file = "3d_raw.out"
      bd = load(file, dir=datapath)
      plotrange = [-50.0, 50.0, -0.5, 0.5]
      X, Z, p = cutdata(bd, "p"; dir="y", sequence=1, plotrange)
      @test p[1] ≈ 0.560976f0
      @test p[2] ≈ 0.53704995f0
      vars = getvars(bd, ["p"])
      @test size(vars["p"]) == (8,8,8)
   end

   @testset "Log" begin
      logfilename = joinpath(datapath, "log_n000001.log")
      head, bd = readlogdata(logfilename)
      @test isa(head, NamedTuple)
      @test extrema(bd) == (-0.105, 258.0)
   end

   @testset "VTK" begin
      file = joinpath(datapath, "3d_bin.dat")
      head, bd, connectivity = readtecdata(file)
      @test maximum(connectivity) ≤ head[:nNode] # check if it's read correctly
      convertTECtoVTU(head, bd, connectivity)
      sha_str = bytes2hex(open(sha1, "out.vtu"))
      @test sha_str == "5b04747666542d802357dec183177f757754a254"
      rm("out.vtu")

      filetag = joinpath(datapath, "3d_mhd_amr/3d__mhd_1_t00000000_n00000000")
      batl = Batl(readhead(filetag*".info"), readtree(filetag)...)
      # local block index check
      @test Batsrus.find_grid_block(batl, [1.0, 0.0, 0.0]) == 4
      @test Batsrus.find_grid_block(batl, [100.0, 0.0, 0.0]) == -100
      connectivity = getConnectivity(batl)
      sha_str = bytes2hex(sha256(string(connectivity)))
      @test sha_str == "c6c5a65a46d86a9ba4096228c1516f89275e45e295cd305eb70c281a770ede74"
   end

   @testset "Plotting" begin
      @testset "Plots" begin
         RecipesBase.is_key_supported(k::Symbol) = true
         # 1D
         file = "1d__raw_2_t25.60000_n00000258.out"
         bd = load(file, dir=datapath)
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bd, "Rho")
         @test getfield(rec[1], 1)[:seriestype] == :path

         file = "z=0_raw_1_t25.60000_n00000258.out"
         bd = load(file, dir=datapath)
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bd, "p")
         @test getfield(rec[1], 1)[:seriestype] == :contourf
      end

      @testset "Makie" begin
         file = "1d__raw_2_t25.60000_n00000258.out"
         bd = load(file, dir=datapath)
         fig, ax, plt = lines(bd, "Rho")
         @test plt isa Lines

         file = "z=0_raw_1_t25.60000_n00000258.out"
         bd = load(file, dir=datapath)
         fig, ax, plt = heatmap(bd, "p")
         @test plt isa Heatmap
      end

      @testset "PyPlot" begin
         # 1D ascii
         file = "1d__raw_2_t25.60000_n00000258.out"
         bd = load(file, dir=datapath, verbose=false)
         plotdata(bd, "p", plotmode="line")
         line = get(gca().lines, 0)
         @test line.get_xdata() ≈ bd.x
         @test line.get_ydata() ≈ bd.w[:,10]

         # 2D structured binary
         file = "z=0_raw_1_t25.60000_n00000258.out"
         bd = load(file, dir=datapath)
         plotdata(bd, "p bx;by", plotmode="contbar streamover")
         @test isa(gca(), PyPlot.PyObject)
         PyPlot.contourf(bd, "p")
         @test isa(gca(), PyPlot.PyObject)
         @test_throws ErrorException PyPlot.contourf(bd, "rho", innermask=true)
         PyPlot.contour(bd, "rho")
         @test isa(gca(), PyPlot.PyObject)
         c = PyPlot.tricontourf(bd, "rho")
         @test c.get_array()[end] == 1.0500000000000003
         PyPlot.tripcolor(bd, "rho")
         @test isa(gca(), PyPlot.PyObject)
         PyPlot.streamplot(bd, "bx;by")
         @test isa(gca(), PyPlot.PyObject)
         p = PyPlot.pcolormesh(bd, "p").get_array()
         @test p[end] == 0.1f0
         plt.close()
         fig = plt.figure()
         ax = fig.add_subplot(111, projection="3d")
         plot_surface(bd, "rho")
         @test isa(gca(), PyPlot.PyObject)
         plt.close()

         # 2D AMR Cartesian
         file = "bx0_mhd_6_t00000100_n00000352.out"
         bd = load(file, dir=datapath)
         plotdata(bd, "P", plotmode="contbar")
         ax = gca()
         @test isa(ax, PyPlot.PyObject)
      end
   end

   @testset "Units" begin
      @test 1.0bu"R" > 2.0bu"Rg"
      file = "y=0_var_1_t00000000_n00000000.out"
      bd = load(file, dir=datapath)
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