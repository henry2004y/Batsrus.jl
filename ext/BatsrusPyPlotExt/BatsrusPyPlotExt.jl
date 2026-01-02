module BatsrusPyPlotExt

using Batsrus
using PyPlot
using Printf
using Interpolations

import Batsrus: meshgrid, findindex, adjust_plotrange!, _resolve_alias
import Batsrus: plotlogdata, cutplot, streamslice, plot_phase

function __init__()
   if VersionNumber(PyPlot.matplotlib.__version__) >= v"3.3"
      PyPlot.matplotlib.rc("image", cmap = "turbo") # set default colormap
   end

   if VersionNumber(PyPlot.matplotlib.__version__) < v"3.5"
      PyPlot.matplotlib.rc("pcolor", shading = "nearest") # newer version default "auto"
   end

   PyPlot.matplotlib.rc("font", size = 16)
   PyPlot.matplotlib.rc("xtick", labelsize = 10)
   PyPlot.matplotlib.rc("ytick", labelsize = 10)
   PyPlot.matplotlib.rc("xtick.minor", visible = true)
   PyPlot.matplotlib.rc("ytick.minor", visible = true)
end

function Batsrus._triangulate_matplotlib(X, Y, W, xi, yi)
   triang = PyPlot.matplotlib.tri.Triangulation(X, Y)
   interpolator = PyPlot.matplotlib.tri.LinearTriInterpolator(triang, W)
   Xi, Yi = meshgrid(xi, yi)
   Wi = interpolator(Xi, Yi) # Always returns Float64!
end

include("pyplot.jl")

end
