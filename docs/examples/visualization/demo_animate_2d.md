# ---
# title: Contour Animation
# id: demo_contour_animation
# date: 2026-04-22
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.12.6
# description: 2D animation using pyplot
# ---

This example shows how to create 2D colored contour animation from series of SWMF outputs using Matplotlib.

The exported `animate` method from the `BatsrusPyPlotExt` extension provides a convenient way to do this. It automatically handles:
* Time-dependent titles.
* Tweakable colorbars.
* Plotting streamlines on top of colored contours by correctly dispatching a streamline into lines and arrows and removing them iteratively frame by frame.

```julia
using Batsrus, PyPlot

filedir = "GM/"

files = filter(file -> startswith(file, "z") && endswith(file, ".out"), readdir(filedir))
# Generate full paths for files
files = joinpath.(filedir, files)

var = "bz"
vmin = -10
vmax = 10
cmap = PyPlot.matplotlib.cm.RdBu_r
streamvars = "bx;by"
density = 1
orientation = "vertical"
outdir = "out/"
extend = "both"

animate(
    files; var, vmin, vmax, outdir, streamvars,
    plot_kwargs = (; cmap),
    stream_kwargs = (; color="white", density),
    cbar_kwargs = (; orientation, extend)
)
```
