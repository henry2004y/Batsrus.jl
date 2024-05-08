# ---
# title: Contour Animation
# id: demo_contour_animation
# date: 2024-05-08
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.3
# description: 2D contour animation using pyplot
# ---

This example shows how to create 2D colored contour animation from series of SWMF outputs.

* Time-dependent title.
* Tweakable colorbar.
* Plotting streamline on top of colored contours is a bit hacky because there is currently no intrinsic methods for removing the streamlines. Our solution here is to dispatch a streamline into lines and arrows and then remove them separately.

```julia
using Batsrus, PyPlot, Printf

"""
    animate2d(files::Vector{String}, var::String; kwargs...)

Saving series of plots of `var` from SWMF output `files`.

# Keywords
- `filedir::String="./"`: input SWMF file directory.
- `outdir::String="out/"`: output directory.
- `vmin=-4`: plot value lower bound.
- `vmax=4`: plot value upper bound.
- `cmap=matplotlib.cm.RdBu_r`: chosen colormap from Matplotlib.
"""
function animate2d(files::Vector{String}, var; filedir="./", outdir="out/", vmin=-4, vmax=4,
   cmap = matplotlib.cm.RdBu_r)
   nfile = length(files)
   x, y = let bd = load(filedir*files[1])
      range(bd.x[1,1,1], bd.x[end,1,1], length=size(bd.x,1)),
      range(bd.x[1,1,2], bd.x[1,end,2], length=size(bd.x,2))
   end
   # colormap norm
   norm = matplotlib.colors.Normalize(vmin, vmax)
   
   fig = plt.figure(figsize=(8, 3), constrained_layout=true)
   ax = plt.axes()

   c = let
      fakedata = zeros(Float32, length(x), length(y))
      ax.pcolormesh(x, y, fakedata'; norm, cmap)
   end
   ax.set_xlabel("x [km]", fontsize=14)
   ax.set_ylabel("y [km]", fontsize=14)
   cb = colorbar(c; ax, orientation="horizontal", location="top", aspect=50)  
   cb.ax.set_ylabel("Bz [nT]", fontsize=14)

   for (i, file) in enumerate(files)
      @info "$i in $nfile"
      bd = load(joinpath(filedir, file))
      c.set_array(bd[var]')
   
      title_str = @sprintf "t = %4.1f s" bd.head.time
      ax.set_title(title_str)
   
      savefig(outdir*lpad(i, 4, '0')*".png", bbox_inches="tight", dpi=200)
   end
   
   close()
end

"""
    animate2d_with_streamline(files::Vector{String}, var::String, streamvars; kwargs...)

Saving series of plots of `var`, together with streamplots of `streamvars`, from SWMF output
`files`.

# Keywords
- `filedir::String="./"`: input SWMF file directory.
- `outdir::String="out/"`: output directory.
- `vmin=-4`: plot value lower bound.
- `vmax=4`: plot value upper bound.
- `cmap=matplotlib.cm.RdBu_r`: chosen colormap from Matplotlib.
"""
function animate2d_with_streamline(files::Vector{String}, var, streamvars="Bx;By";
   filedir="./", outdir="out/", vmin=-4, vmax=4, cmap=matplotlib.cm.RdBu_r)
   nfile = length(files)
   x, y = let bd = load(filedir*files[1])
      range(bd.x[1,1,1], bd.x[end,1,1], length=size(bd.x,1)),
      range(bd.x[1,1,2], bd.x[1,end,2], length=size(bd.x,2))
   end
   # colormap norm
   norm = matplotlib.colors.Normalize(vmin, vmax)

   fig = plt.figure(figsize=(8, 3), constrained_layout=true)
   ax = plt.axes()

   @info "1 in $nfile"
   c = ax.pcolormesh(x, y, bd[var]'; norm, cmap)

   ax.set_xlabel("x [km]", fontsize=14)
   ax.set_ylabel("y [km]", fontsize=14)
   cb = colorbar(c; ax, orientation="horizontal", location="top", aspect=50)  
   cb.ax.set_ylabel("Bz [nT]", fontsize=14)
   title_str = @sprintf "t = %4.1f s" bd.head.time
   ax.set_title(title_str)

   st = streamplot(bd, streamvars, ax; color="gray", density=2)

   save_and_clean!(1, outdir, st, ax)

   for i in 2:nfile
      @info "$i out of $nfile"
      bd = load(joinpath(filedir, files[i]))
      c.set_array(bd["Bz"]')

      st = streamplot(bd, streamvars, ax; color="gray", density=2)
      ax.set_xlim(x[1], x[end])
      ax.set_ylim(y[1], y[end])

      title_str = @sprintf "t = %4.1f s" bd.head.time
      ax.set_title(title_str)
   
      save_and_clean!(i, outdir, st, ax)
   end
   
   close()
end

function save_and_clean!(i::Int, outdir::String, st, ax)
   savefig(outdir*lpad(i, 4, '0')*".png", bbox_inches="tight", dpi=200)
   # Clean up streamlines
   st.lines.remove()
   for art in ax.get_children()
      if !pybuiltin(:isinstance)(art, matplotlib.patches.FancyArrowPatch)
         continue
      end
      art.remove()
   end
end

####################
# Data directory
filedir = "./"
# Plot variables
var = "By"

# Find simulation data
files = let
   pick_file = file -> startswith(file, "z") && endswith(file, ".out")
   filter(pick_file, readdir(filedir))
end

animate2d(files, var)
#animate2d_with_streamline(files, var, "Bx;By")
```
