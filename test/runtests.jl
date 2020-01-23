# Test of VisAna

ENV["MPLBACKEND"]="agg" # no GUI

using VisAna, PyPlot, Test

@info("VisAna is calling PyPlot for plotting.")

@testset "reading 1D ascii" begin
   filename = "1d__raw_2_t25.60000_n00000258.out"
   head, data, list = readdata(filename, verbose=true)
   @test isa(head[1], Dict)
   @test isa(data[1], Data)
   @test isa(list[1], FileList)
   plotdata(data[1],head[1],"p",plotmode="line")
   line = gca().lines[1]
   @test line.get_xdata() ≈ data[1].x
   @test line.get_ydata() ≈ data[1].w[:,10]
end

@testset "reading 2D structured binary" begin
   filename = "z=0_raw_1_t25.60000_n00000258.out"
   head, data, list = readdata(filename)
   @test isa(head[1], Dict)
   @test isa(data[1], Data)
   @test isa(list[1], FileList)
   plotdata(data[1],head[1],"p bx;by",plotmode="contbar streamover")
   ax = gca()
   @test isa(ax, PyPlot.PyObject)
   contourf(data[1],head[1],"p")
   ax = gca()
   @test isa(ax, PyPlot.PyObject)
end

@testset "reading 2D unstructured binary" begin
   #filename = "z=0_raw_1_t25.60000_n00000258.out"
   #head, data, list = readdata(filename)
end

@testset "reading 3D structured binary" begin
   #filename = "3d_structured.out";
   #head, data, list = readdata(filename,verbose=false);
end

@testset "log" begin
   logfilename = "log_n000001.log"
   head, data = readlogdata(logfilename)
   @test isa(head, Dict)
   @test isa(data, Array)
end

@testset "vtk" begin
   @info("VTK conversion test.")
   filename = "3d_bin.dat"
   head, data, connectivity  = readtecdata(filename, IsBinary=true)
   @test maximum(connectivity) ≤ head[:nNode] # check if it's read correctly
   convertVTK(head, data, connectivity)
   @test isfile("3DBATSRUS.vtu")
   rm("3DBATSRUS.vtu")
end

@testset "MVA" begin
   @info("Minimum variance analysis test.")
   include("../src/space/MVA.jl")
   F = MVA_analysis("Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
   eigenRef = [2079.360,  78.142,   34.309]
   @test F.values ≈ eigenRef atol=1e-3
end
