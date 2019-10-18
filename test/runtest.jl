# Test of VisAna
#
# Hongyang Zhou, hyzhou@umich.edu 07/30/2019

# For better usage in REPL, I need to provide a script or wrapper to the full
# functions. For example, I want something like
# readdata
# in IDL. To achieve that, probably I need a script to do so.

using VisAna, Test

#@testset "some tests" begin
#   @test filehead, data, filelist = readdata(filename);
#end

# 1D
filename = "1d_bin.out";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"p",plotmode="line")
plotdata(data[1],filehead[1],"p",plotmode="linegrid")
animatedata(data[1],filehead[1],"p")

# 2D Cartesian (structured)
filename = "z=0_raw_1_t10.07620_n00000102.out";
filehead, data, filelist = readdata(filename,verbose=false);

plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar streamover")
plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar quiverover")
plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar streamover",
   density=2.0)
plotdata(data[1],filehead[1],"p",plotmode="grid")
plotdata(data[1],filehead[1],"p",plotmode="contbar",
   plotrange=[-50., 50., -1., 1.])
plotdata(data[1],filehead[1],"p",plotmode="contbar")
plotdata(data[1],filehead[1],"p",plotmode="contbarlog")
plotdata(data[1],filehead[1],"p",plotmode="surfbar")

# 2D unstructured
filename = "y=0_unstructured.outs";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar",
   plotrange=[-2.2, 2.2, -2.2, 2.2])
plotdata(data[1], filehead[1], "rho", plotmode="trimesh")
plotdata(data[1], filehead[1], "rho", plotmode="tricont")

animatedata(data[1],filehead[1],filelist[1],"rho",plotmode="contbar")

# 2D structured spherical ???
filename = "y_structured.out"; #???
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar")

# Similar API as matplotlib

# 2D contour
contour(data[1], filehead[1], "p")

# 2D contourf
contourf(data[1], filehead[1], "p")
contourf(data[1], filehead[1], "p", levels, plotrange=[-10,10,-Inf,Inf],
   plotinterval=0.1)

# surface
plot_surface(data[1], filehead[1], "p")

# triangle filled contour
tricontourf(data[1], filehead[1], "p")

# triangle surface plot
plot_trisurf(data[1], filehead[1], "p")

# 3D box
filename = "3d_box.out";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1], filehead[1], "ps1 bx;bz", plotmode="contbar stream", cut="y",
   cutPlaneIndex=65)
plotdata(data[1], filehead[1], "ps1 bx;by", plotmode="contbar stream", cut="z")
plotdata(data[1], filehead[1], "bx", plotmode="contbar", cut="y")
plotdata(data[1], filehead[1], "bx", plotmode="contbar", cut="y",
   plotrange=[-1.4,-1.1,0.70,0.78])

# 3D structured spherical
filename = "3d_structured.out";
filehead, data, filelist = readdata(filename,verbose=false);

# log file
logfilename = "shocktube.log";
filehead, data = readlogdata(logfilename)

# Tecplot ascii to vtk conversion
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
tec = readtecdata.(filenames, false) # head, data, connectivity
for (i, outname) in enumerate(filenames)
   convertVTK(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])
end



using PyCall, PyPlot

include("VisAna.jl")
using .VisAna

filename = "y*.outs"
filehead, data, filelist = readdata(filename,npict=20,verbose=false);

X = vec(data[1].x[:,:,1])
Y = vec(data[1].x[:,:,2])
W = vec(data[1].w[:,:,10])

# Perform linear interpolation of the data (x,y) on grid(xi,yi)
plotrange = zeros(4)
plotrange[1] = -3 # minimum(X)
plotrange[2] = 3  # maximum(X)
plotrange[3] = -3 # minimum(Y)
plotrange[4] = 3  # maximum(Y)

#plotinterval = X[2] - X[1]
plotinterval = 0.01

xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

triang = matplotlib.tri.Triangulation(X,Y)
interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
np = pyimport("numpy")
Xi, Yi = np.meshgrid(xi, yi)
wi = interpolator(Xi, Yi)

# mask a circle in the middle:
interior = (Xi.^2 .+ Yi.^2 .< 1.0);
wi[interior] .= np.ma.masked

c = contourf(xi,yi,wi,50)
colorbar()
using MATLAB
px, py = mxcall(:gradient, 2, wi)
