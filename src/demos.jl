# Test of VisAna
#
# Hongyang Zhou, hyzhou@umich.edu 07/30/2019

# For better usage in REPL, I need to provide a script or wrapper to the full
# functions. For example, I want something like
# readdata
# in IDL. To achieve that, probably I need a script to do so.

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

####################
using VisAna, PyPlot

include("../src/trace2d.jl")

filename = "y=0_var_1_t00000000_n00000000.out"
head, data, list = readdata(filename,dir="test")

#streamplot(data[1], head[1], "bx;bz")

bx = data[1].w[:,:,5]
bz = data[1].w[:,:,7]
x  = data[1].x[:,1,1]
z  = data[1].x[1,:,2]

#xstart = -120.0:10.0:10.
#zstart = fill(10.0, size(xstart))

seeds = select_seeds(x,z)

for i = 1:size(seeds)[2]
   xs = seeds[1,i]
   zs = seeds[2,i]
   # forward
   x1, y1 = trace2d_eul(bx, bz, xs, zs, x, z, ds=0.1, maxstep=1000, gridType="ndgrid")
   plot(x1,y1,"--")
   # backward
   x2, y2 = trace2d_rk4(-bx, -bz, xs, zs, x, z, ds=0.1, maxstep=1000, gridType="ndgrid")
   plot(x2,y2,"-")
end
axis("equal")
