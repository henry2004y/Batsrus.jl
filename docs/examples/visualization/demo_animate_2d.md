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
* Consistent colors across frames: if `vmin`/`vmax` are not given, **all** snapshots are scanned for a global color range (frames are cached during the scan, so each file is only read once).
* AMR-aware sampling: with the default `plotinterval = nothing`, the finest grid resolution is auto-detected from the first snapshot (see [`finest_resolution`](@ref)), which avoids oversampling or undersampling AMR and unstructured output.
* Automatic colormap: bipolar data gets a symmetric `RdBu_r` scale, non-negative data gets `turbo`. Override with `colormap` or `plot_kwargs`.

```julia
using Batsrus, PyPlot

filedir = "GM/"

files = filter(file -> startswith(file, "z") && endswith(file, ".out"), readdir(filedir))
# Generate full paths for files
files = joinpath.(filedir, files)

var = "bz"
streamvars = "bx;by"
outdir = "out/"

# vmin/vmax and the colormap are determined automatically; bipolar data such as
# bz gets a symmetric RdBu_r color scale consistent across all frames.
animate(
    files; var, outdir, streamvars,
    stream_kwargs = (; color = "white", density = 1),
    cbar_kwargs = (; orientation = "vertical", extend = "both")
)
```
