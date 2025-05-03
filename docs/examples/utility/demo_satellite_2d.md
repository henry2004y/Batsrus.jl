# ---
# title: Virtual Satellite Extraction
# id: demo_satellite_2d
# date: 2025-05-03
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.11.4
# description: Satellite extraction from 2D slices
# ---

This example shows how to extract virtual satellite data from FLEKS output as a time series.

```julia
using Batsrus, JLD2

function save_static_satellite_data(
   loc;
   filedir = "./",
   outname = "satellite.jld2",
   tstep = 1.0,
)
   pick_file = file -> startswith(file, "z") && endswith(file, ".out")

   files = [joinpath(filedir, f) for f in readdir(filedir) if pick_file(f)]

   trange, v = get_timeseries(files, loc; tstep)
   jldsave(outname; trange, v, loc)
end

#####
mi2me = 16
f = 2
deg = 0
dim = 3
dx = 1 # [de]
tstep = 1.0 # output time cadence
# Data directory
topdir = "/home1/"
filedir = joinpath(topdir, "RESULT/run1/PC/")
outdir = joinpath(topdir, "run1/")

if !isdir(outdir)
   println("Creating figure saving directory:")
   println(outdir)
   mkdir(outdir)
else
   println("Output data saving directory:")
   println(outdir)
end

locs = [Float32[i, 0] for i in range(10, 20, step=5)]

for i in eachindex(locs)
   @info "Location $(i)"
   loc = locs[i]
   outname = joinpath(outdir, "satellite_$(dim)d$(deg)deg_mi2me$(mi2me)_f$(f)_$i.jld2")
   isfile(outname) && continue
   save_static_satellite_data(loc; filedir, outname, tstep)
end
```
