# Test of VisAna
#
# Hongyang Zhou, hyzhou@umich.edu 07/30/2019

# For better usage in REPL, I need to provide a script or wrapper to the full
# functions. For example, I want something like
# readdata
# in IDL. To achieve that, probably I need a script to do so.

using Test

include("VisAna.jl")
using .VisAna

filename = "1d_bin.out";

@testset "some tests" begin
   @test filehead, data, filelist = readdata(filename);
end

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
plotdata(data[1],filehead[1],"rho",plotmode="contbar")
plotdata(data[1],filehead[1],"rho",plotmode="trimesh")
plotdata(data[1],filehead[1],"rho",plotmode="tricont")

animatedata(data[1],filehead[1],filelist[1],"rho",plotmode="contbar")

# 2D structured spherical ???
filename = "y_structured.out"; #???
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar")


# 3D box
filename = "box.outs";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"bx",plotmode="contbar",cut="y")
plotdata(data[1],filehead[1],"bx",plotmode="contbar",cut="y",
   plotrange=[-1.4,-1.1,0.70,0.78])

# 3D structured spherical
filename = "3d_structured.out";
filehead, data, filelist = readdata(filename,verbose=false);

# log file
logfilename = "shocktube.log";
filehead, data = readlogdata(logfilename)


using PyCall, PyPlot

include("VisAna.jl")
using .VisAna

filename = "y*.outs"
filehead, data, filelist = readdata(filename,verbose=false);

X = vec(data[1].x[:,:,1])
Y = vec(data[1].x[:,:,2])
W = vec(data[1].w[:,:,8])

# Perform linear interpolation of the data (x,y) on grid(xi,yi)
plotrange = zeros(4)
plotrange[1] = minimum(X)
plotrange[2] = maximum(X)
plotrange[3] = minimum(Y)
plotrange[4] = maximum(Y)

#plotinterval = X[2] - X[1]
plotinterval = 1.0

xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

triang = matplotlib.tri.Triangulation(X,Y)
interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
np = pyimport("numpy")
Xi, Yi = np.meshgrid(xi, yi)
wi = interpolator(Xi, Yi)

c = contourf(xi,yi,wi)
using MATLAB
px, py = mxcall(:gradient, 2, wi)


include("VisAna.jl")
using .VisAna
filename = "3d.dat"
filehead, data, filelist = readdata(filename,verbose=false);

X = vec(data[1,:]);
Y = vec(data[2,:]);
Z = vec(data[3,:]);
W = vec(data[8,:]);

using PyCall

np = pyimport("numpy")
scipy = pyimport("scipy")
interpolate = pyimport("scipy.interpolate")
griddata = interpolate.griddata

s = pybuiltin(:slice)
XYZ = get(np.mgrid, (s(0,1,10im), s(0,1,20im), s(0,1,10im)))
grid_x, grid_y, grid_z = XYZ[1,:,:,:], XYZ[2,:,:,:], XYZ[3,:,:,:]

points = np.random.rand(1000, 3)
values = py"func"(points[:,1], points[:,2], points[:,3])

grid_z0 = griddata(points, values, (grid_x, grid_y, grid_z), method="nearest")
grid_z1 = griddata(points, values, (grid_x, grid_y, grid_z), method="linear")
grid_z2 = griddata(points, values, (grid_x, grid_y, grid_z), method="cubic")
