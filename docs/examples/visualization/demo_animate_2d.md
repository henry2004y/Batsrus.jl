# ---
# title: Contour Animation
# id: demo_contour_animation
# date: 2024-05-08
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.3
# description: 2D contour animation using pyplot
# ---

This example shows how to create 2D colored contour animation from series of SWMF outputs using Matplotlib.

* Time-dependent title.
* Tweakable colorbar.
* Plotting streamline on top of colored contours is a bit hacky because there is currently no intrinsic methods for removing the streamlines. Our solution here is to dispatch a streamline into lines and arrows and then remove them separately.

```julia
using Batsrus, PyPlot, Printf, PyCall

const fontsize = 14

"""
    animate(files::Vector{String}; kwargs...)

Generate figures of colored contour, optionally with streamlines.

# Keywords
- `var`: variable to plot with pcolormesh.
- `vmin`: minimum plotting value.
- `vmax`: maximum plotting value.
- `plotrange`: 2D plotting spatial range.
- `plotinterval`: spatial sampling interval.
- `cmap`: colormap accepted by Matplotlib.
- `streamvars`: if set, add streamlines to the colored contour.
- `density`: streamline density.
- `orientation`: colorbar orientation, "horizontal" or "vertical".
- `extend`: colorbar extension, "neither", "both", "min", "max".
- `use_conventional_notation`: if true, denote the vertical axis as "z".
- `broken_streamlines`: if true, streamlines can be broken.
- `overwrite`: if true, overwrite the existing image files.
- `filedir`: source data directory.
- `outdir`: output image directory.
"""
function animate(files::Vector{String}; var="Uy", vmin=-Inf, vmax=Inf,
   plotrange=[-Inf, Inf, -Inf, Inf], plotinterval=Inf,
   cmap=matplotlib.cm.RdBu_r, streamvars=nothing, density=1, orientation="horizontal",
   extend="neither", use_conventional_notation=false, broken_streamlines=false,
   overwrite=false, filedir="./", outdir="out/")

   bd = load(filedir*files[1])

   nfile = length(files)
   if isnan(vmin+vmax)
      vmin = isinf(vmin) ? minimum(v) : vmin
      vmax = isinf(vmax) ? maximum(v) : vmax
   end

   norm = matplotlib.colors.Normalize(vmin, vmax)

   filetype =
      if occursin("region", files[1])
         :PC
      else
         :GM
      end

   if use_conventional_notation # switch y and z axis
      if var[2] == 'y'
         cb_str = 
            if startswith(var, "U")
               "Uz [km/s]"
            elseif startswith(var, "B")
               "Bz [nT]"
            elseif startswith(var, "E")
               L"Ez [$\mu\mathrm{V}/\mathrm{m}$]"
            else
               var
            end
      elseif var[2] == 'z'
         cb_str = 
            if startswith(var, "U")
               "Uy [km/s]"
            elseif startswith(var, "B")
               "By [nT]"
            elseif startswith(var, "E")
               L"Ey [$\mu\mathrm{V}/\mathrm{m}$]"
            else
               var
            end
      end
   else
      cb_str = 
         if startswith(var, "U")
            var * " [km/s]"
         elseif startswith(var, "B")
             var * " [nT]"
         elseif startswith(var, "E")
            var * L" [\mu\mathrm{V}/\mathrm{m}]"
         else
            var
         end
   end

   fig = plt.figure(figsize=(4, 8), constrained_layout=true)
   ax = plt.axes()
   ax.set_xlabel(L"x [$\mathrm{R}_\mathrm{E}$]"; fontsize)
   ax.set_ylabel(L"y [$\mathrm{R}_\mathrm{E}$]"; fontsize)

   @info "1 out of $nfile"
   x, y, w = Batsrus.getdata2d(bd, var, plotrange, plotinterval; innermask=false)
   c = ax.pcolormesh(x, y, w; norm, cmap)
   ax.set_aspect("equal", adjustable="box", anchor="C")

   if orientation == "horizontal"
      cb = colorbar(c; ax, orientation, location="top", aspect=50, extend)
   else
      cb = colorbar(c; ax, orientation, extend)  
   end
   cb.ax.set_ylabel(cb_str; fontsize)
   #cb.ax.set_title(cb_str; fontsize)
   title_str = @sprintf "t = %4.1f s" bd.head.time
   ax.set_title(title_str)

   if bd.head.gencoord
      if !isnothing(streamvars)
         xi, yi, v1, v2 = get_vector(bd, streamvars, plotrange, plotinterval)
         st = ax.streamplot(xi, yi, v1, v2; color="gray", density)
      end
   else
      if !isnothing(streamvars)
         st = streamplot(bd, streamvars, ax; color="gray", density, broken_streamlines)
      end
   end

   savefig(outdir*lpad(1, 4, '0')*".png", bbox_inches="tight", dpi=200)

   if !isnothing(streamvars) # Clean up streamlines
      st.lines.remove()
      for art in ax.get_children()
         if !pybuiltin(:isinstance)(art, matplotlib.patches.FancyArrowPatch)
            continue
         end
         art.remove()
      end
   end

   for i in 2:nfile
      @info "$i out of $nfile"
      outname = outdir*lpad(i, 4, '0')*".png"
      if !overwrite
         isfile(outname) && continue
      end
      bd = load(joinpath(filedir, files[i]))
      if bd.head.gencoord
         _, _, wi = getdata2d(bd, var, plotrange, plotinterval; innermask=false)
      else
         wi = bd[var]'
      end

      c.set_array(wi)

      if !isnothing(streamvars)
         if bd.head.gencoord
            xi, yi, v1, v2 = get_vector(bd, streamvars, plotrange, plotinterval)
            st = ax.streamplot(xi, yi, v1, v2; color="gray", density)
         else
            st = streamplot(bd, "Bx;By", ax; color="gray", density, broken_streamlines)
         end
      end

      title_str = @sprintf "t = %4.1f s" bd.head.time
      ax.set_title(title_str)
   
      savefig(outname, bbox_inches="tight", dpi=200)
      if !isnothing(streamvars) # Clean up streamlines
         st.lines.remove()
         for art in ax.get_children()
            if !pybuiltin(:isinstance)(art, matplotlib.patches.FancyArrowPatch)
               continue
            end
            art.remove()
         end
      end
   end
   
   close()
end

function get_vector(bd, var, plotrange, plotinterval)
   x, w = bd.x, bd.w
   X, Y = eachslice(x, dims=3)
   X, Y = vec(X), vec(Y)
   varstream = split(var, ";")
   var1_ = findfirst(x->lowercase(x)==lowercase(varstream[1]), bd.head.wnames)
   var2_ = findfirst(x->lowercase(x)==lowercase(varstream[2]), bd.head.wnames)
   # Create grid values first.
   xi = range(Float64(plotrange[1]), stop=Float64(plotrange[2]), step=plotinterval)
   yi = range(Float64(plotrange[3]), stop=Float64(plotrange[4]), step=plotinterval)

   tr = matplotlib.tri.Triangulation(X, Y)
   Xi, Yi = Batsrus.meshgrid(xi, yi)

   W1 = @views w[:,:,var1_] |> vec
   interpolator = matplotlib.tri.LinearTriInterpolator(tr, W1)
   v1 = interpolator(Xi, Yi)

   W2 = @views w[:,:,var2_] |> vec
   interpolator = matplotlib.tri.LinearTriInterpolator(tr, W2)
   v2 = interpolator(Xi, Yi)

   xi, yi, v1, v2
end

#######################
filedir = "GM/"

files = filter(file -> startswith(file, "z") && endswith(file, ".out"), readdir(filedir))

var = "Bz"
vmin = -10
vmax = 10
cmap = matplotlib.cm.RdBu_r
streamvars = "Bx;By"
density = 1
orientation = "vertical"
outdir = "out/"
extend = "both"
use_conventional_notation = true

animate(files; filedir, outdir, var, vmin, vmax, orientation, cmap, streamvars, density,
   extend, use_conventional_notation)
```
