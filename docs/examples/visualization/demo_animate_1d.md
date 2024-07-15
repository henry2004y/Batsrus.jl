# ---
# title: Line Animation
# id: demo_line_animation
# date: 2024-05-08
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.3
# description: 1D line animation using pyplot
# ---

This example shows how to create 1D space line animation from series of SWMF outputs.

* Vertical line moving at the solar wind speed for reference.
* Twin axes sharing the same x-axis for showing two quantities in one frame.
* Time-dependent title.

```julia
using Batsrus, PyPlot, Printf

"""
    animate1d(files::Vector{String}, vars::Vector{String}; kwargs...)

Saving series of plots of `vars` from SWMF output `files`.

# Keywords
- `filedir::String="./"`: input SWMF file directory.
- `outdir::String="out/"`: output directory.
- `output_interval::Int=1`: Timestep interval for output files.
- `vmin=-1`: plot value lower bound.
- `vmax=1`: plot value upper bound.
- `refloc=1e5`: reference vertical line initial location.
"""
function animate1d(files::Vector{String}, vars::Vector{String};
   filedir::String="./", outdir::String="out/", output_interval::Int=1,
   vmin=-1, vmax=1, refloc=1e5)
   nfile = length(files)
   x = let
      bd = load(joinpath(filedir, files[1]))
      bd.x[:,1,1]
   end

   fig = plt.figure(figsize=(12,5), constrained_layout=true)
   ax = plt.axes(xlim=extrema(x), ylim=(vmin, vmax))

   line, = ax.plot([], [], lw=1)
   # Move reference line forward in the canvas
   vl = ax.axvline(0, ls="-", color="tab:brown", lw=1, zorder=10)

   color = "tab:blue"
   ax.set_xlabel("x [km]"; fontsize=14)
   ax.set_ylabel(var*" [nT]"; fontsize=14, color)
   ax.tick_params(axis="y", labelcolor=color)

   ax2 = ax.twinx()
   ax2.set_ylim(-600, 0)

   color = "tab:red"
   ax2.set_ylabel("Ux [km/s]"; fontsize=14, color)

   line2, = ax2.plot([], []; color, alpha=0.7)
   ax2.tick_params(axis="y", labelcolor=color)


   for (i, file) in enumerate(files)
      @info "$i in $nfile"
      outname = outdir*lpad(i, 4, '0')*".png"
      isfile(outname) && continue

      bd = load(joinpath(filedir, file))

      d = bd[vars[1]][:,1]
      title_str = @sprintf "t = %4.1f s" bd.head.time
      line.set_data(x, d)
      d2 = bd[vars[2]][:,1]
      line2.set_data(x, d2)

      ax.set_title(title_str)

      refloc -= output_interval * VSW
      if refloc <= 0
         refloc += 1e5
      end
      vl.set_xdata([refloc, refloc])

      savefig(outname, bbox_inches="tight", dpi=200)
   end

   close()
end

####################
# Constants
const VSW = 500.0  # Solar wind speed, [km/s]

# Data directory
filedir = "./"
# Plot variables
vars = ["By", "uxs1"]

# Find simulation data
files = let
   pick_file = file -> startswith(file, "z") && endswith(file, ".out")
   filter(pick_file, readdir(filedir))
end

animate1d(files, vars)
```
