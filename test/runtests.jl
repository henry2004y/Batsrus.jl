# Tests of BATSRUS.jl

using Batsrus, Test, SHA, LazyArtifacts
using Batsrus.UnitfulBatsrus, Unitful

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
      data = load(file, dir=datapath, verbose=true)
      @test startswith(repr(data), "filename : 1d")
      @test Batsrus.setunits(data.head, "NORMALIZED")
      @test isa(data.head, NamedTuple)
      @test extrema(data.x) == (-127.5, 127.5)
      @test extrema(data.w) == (-0.79960780498, 1.9394335293)
   end

   @testset "Reading 2D structured binary" begin
      file = "z=0_raw_1_t25.60000_n00000258.out"
      data = load(file, dir=datapath)
      @test data.head.time == 25.6f0
      @test extrema(data.x) == (-127.5f0, 127.5f0)
      @test extrema(data.w) == (-0.79985905f0, 1.9399388f0)
   end

   @testset "Reading 2D unstructured binary" begin
      #file = "z=0_raw_1_t25.60000_n00000258.out"
      #data = load(file)
   end

   @testset "Reading 3D structured binary" begin
      file = "3d_raw.out"
      data = load(file, dir=datapath)
      plotrange = [-50.0, 50.0, -0.5, 0.5]
      X, Z, p = cutdata(data, "p"; dir="y", sequence=1, plotrange)
      @test p[1] ≈ 0.560976f0
      @test p[2] ≈ 0.53704995f0
      vars = getvars(data, ["p"]) 
      @test size(vars["p"]) == (8,8,8) 
   end

   @testset "Log" begin
      logfilename = joinpath(datapath, "log_n000001.log")
      head, data = readlogdata(logfilename)
      @test isa(head, NamedTuple)
      @test extrema(data) == (-0.105, 258.0)
   end

   @testset "VTK" begin
      file = joinpath(datapath, "3d_bin.dat")
      head, data, connectivity = readtecdata(file)
      @test maximum(connectivity) ≤ head[:nNode] # check if it's read correctly
      convertTECtoVTU(head, data, connectivity)
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
      using PyPlot
      ENV["MPLBACKEND"]="agg" # no GUI
      @testset "1D ascii" begin
         file = "1d__raw_2_t25.60000_n00000258.out"
         data = load(file, dir=datapath, verbose=false)
         plotdata(data, "p", plotmode="line")
         line = get(gca().lines, 0)
         @test line.get_xdata() ≈ data.x
         @test line.get_ydata() ≈ data.w[:,10]
      end

      @testset "2D structured binary" begin
         file = "z=0_raw_1_t25.60000_n00000258.out"
         data = load(file, dir=datapath)
         plotdata(data, "p bx;by", plotmode="contbar streamover")
         @test isa(gca(), PyPlot.PyObject)
         contourf(data, "p")
         @test isa(gca(), PyPlot.PyObject)
         @test_throws ErrorException contourf(data, "rho", innermask=true)
         contour(data, "rho")
         @test isa(gca(), PyPlot.PyObject)
         tricontourf(data, "rho")
         @test isa(gca(), PyPlot.PyObject)
         streamplot(data, "bx;by")
         @test isa(gca(), PyPlot.PyObject)
         plt.close()
         fig = plt.figure()
         ax = fig.add_subplot(111, projection="3d")
         plot_surface(data, "rho")
         @test isa(gca(), PyPlot.PyObject)
         plt.close()
      end

      @testset "2D AMR Cartesian" begin
         file = "bx0_mhd_6_t00000100_n00000352.out"
         data = load(file, dir=datapath)
         plotdata(data, "P", plotmode="contbar")
         ax = gca()
         @test isa(ax, PyPlot.PyObject)
      end
   end

   @testset "Units" begin
      @test 1.0bu"R" > 2.0bu"Rg"
      file = "y=0_var_1_t00000000_n00000000.out"
      data = load(file, dir=datapath)
      varunit = getunit(data, "Rho")
      @test varunit == bu"amucc"
      varunit = getunit(data, "Ux")
      @test varunit == u"km/s"
      varunit = getunit(data, "Bx")
      @test varunit == u"nT"
      varunit = getunit(data, "P")
      @test varunit == u"nPa"
      varunit = getunit(data, "jx")
      @test dimension(varunit) == dimension(Unitful.A/Unitful.m^2)
      varunit = getunit(data, "ex")
      @test varunit == u"mV/m"
   end

end