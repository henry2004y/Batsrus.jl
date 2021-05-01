# Tests of BATSRUS.jl

using Batsrus, Test, SHA
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
   @testset "Reading 1D ascii" begin
      filename = "1d__raw_2_t25.60000_n00000258.out"
      data = readdata(filename, dir="data", verbose=true)
      @test isa(data.head, NamedTuple)
      @test extrema(data.x) == (-127.5, 127.5)
      @test extrema(data.w) == (-0.79960780498, 1.9394335293)
   end

   @testset "Reading 2D structured binary" begin
      filename = "z=0_raw_1_t25.60000_n00000258.out"
      data = readdata(filename, dir="data")
      @test data.head.time == 25.6f0
      @test extrema(data.x) == (-127.5f0, 127.5f0)
      @test extrema(data.w) == (-0.79985905f0, 1.9399388f0)
   end

   @testset "Reading 2D unstructured binary" begin
      #filename = "z=0_raw_1_t25.60000_n00000258.out"
      #data = readdata(filename)
   end

   @testset "Reading 3D structured binary" begin
      filename = "3d_raw.out"
      data = readdata(filename, dir="data")
      plotrange = [-50.0, 50.0, -0.5, 0.5]
      X, Z, p = cutdata(data, "p", cut='y', cutPlaneIndex=1, plotrange=plotrange)
      @test p[1] ≈ 0.560976f0
      @test p[2] ≈ 0.53704995f0
      vars = getvars(data, ["p"]) 
      @test size(vars["p"]) == (8,8,8) 
   end

   @testset "Log" begin
      logfilename = "data/log_n000001.log"
      head, data = readlogdata(logfilename)
      @test isa(head, NamedTuple)
      @test extrema(data) == (-0.105, 258.0)
   end

   @testset "VTK" begin
      filename = "data/3d_bin.dat"
      head, data, connectivity = readtecdata(filename)
      @test maximum(connectivity) ≤ head[:nNode] # check if it's read correctly
      convertTECtoVTU(head, data, connectivity)
      sha_str = bytes2hex(open(sha1, "out.vtu"))
      @test sha_str == "5b04747666542d802357dec183177f757754a254"
      rm("out.vtu")

      filetag = "data/3d_mhd_amr/3d__mhd_1_t00000000_n00000000"
      run(`tar -C data -zxf data/3d_mhd_amr.tar.gz`)
      batl = Batl(readhead(filetag*".info"), readtree(filetag)...)
      # local block index check
      @test Batsrus.find_grid_block(batl, [1.0, 0.0, 0.0]) == 4 
      @test Batsrus.find_grid_block(batl, [100.0, 0.0, 0.0]) == -100
      connectivity = getConnectivity(batl)
      sha_str = bytes2hex(sha256(string(connectivity)))
      @test sha_str == "c6c5a65a46d86a9ba4096228c1516f89275e45e295cd305eb70c281a770ede74"
      rm("data/3d_mhd_amr", recursive=true)
   end

   @testset "Plotting" begin
      using PyPlot
      ENV["MPLBACKEND"]="agg" # no GUI
      @testset "Plotting 1D ascii" begin
         filename = "1d__raw_2_t25.60000_n00000258.out"
         data = readdata(filename, dir="data", verbose=false)
         plotdata(data, "p", plotmode="line")
         line = gca().lines[1]
         @test line.get_xdata() ≈ data.x
         @test line.get_ydata() ≈ data.w[:,10]
      end

      @testset "Plotting 2D structured binary" begin
         filename = "z=0_raw_1_t25.60000_n00000258.out"
         data = readdata(filename, dir="data")
         plotdata(data,"p bx;by", plotmode="contbar streamover")
         ax = gca()
         @test isa(ax, PyPlot.PyObject)
         contourf(data,"p")
         ax = gca()
         @test isa(ax, PyPlot.PyObject)
      end
   end

   @testset "Units" begin
      @test 1.0bu"R" > 2.0bu"Rg"
      filename = "y=0_var_1_t00000000_n00000000.out"
      data = readdata(filename, dir="data")
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