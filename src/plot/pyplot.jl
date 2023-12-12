# Plotting functionalities.

using PyPlot
using Interpolations: cubic_spline_interpolation

export plotdata, plotlogdata, plot, scatter, contour, contourf, plot_surface, tripcolor,
   tricontourf, plot_trisurf, streamplot, streamslice, quiver, cutplot, pcolormesh

@static if matplotlib.__version__ >= "3.3"
   matplotlib.rc("image", cmap="turbo") # set default colormap
end

@static if matplotlib.__version__ < "3.5"
   matplotlib.rc("pcolor", shading="nearest") # newer version default "auto"
end

matplotlib.rc("font", size=14)
matplotlib.rc("xtick", labelsize=10)
matplotlib.rc("ytick", labelsize=10)
matplotlib.rc("xtick.minor", visible=true)
matplotlib.rc("ytick.minor", visible=true)

"""
    plotlogdata(data, head, func; plotmode="line")

Plot information from log file.
# Input arguments
- `data::Array`: output data.
- `head::NamedTuple`: header info.
- `func::String`: variables for plotting.
- `plotmode::String`: type of plotting ["line","scatter"].
"""
function plotlogdata(data, head::NamedTuple, func::AbstractString; plotmode="line")
   vars     = split(func)
   plotmode = split(plotmode)

   for (ivar, var) in enumerate(vars)
      varIndex_ = findfirst(x->x==lowercase(var), lowercase.(head.variables))
      isnothing(varIndex_) && error("$(var) not found in file header variables!")

      figure()
      if plotmode[ivar] == "line"
         plot(data[1,:],data[varIndex_,:])
      elseif plotmode[ivar] == "scatter"
         scatter(data[1,:], data[varIndex_,:])
      else
         throw(ArgumentError("unknown plot mode $(plotmode[ivar])!"))
      end
      xlabel(head.variables[1])
      ylabel(head.variables[varIndex_])
      title("log file data")
   end

end


"""
    plotdata(data, func, args, kwargs...)

Plot the variable from SWMF output.

`plotdata(data, "p", plotmode="contbar")`

`plotdata(data, "p", plotmode="grid")`

`plotdata(data, func, plotmode="trimesh", plotrange=[-1.0, 1.0, -1.0, 1.0], plotinterval=0.2)`

# Arguments
- `bd::BATLData`: BATSRUS data to be visualized.
- `func::String`: variables for plotting.

# Keywords
- `plotmode::String`: type of plotting ["cont","contbar"]...
- `plotrange::Vector`: range of plotting.
- `plotinterval`: interval for interpolation.
- `levels`: levels of contour.
- `innermask`: Bool for masking a circle at the inner boundary.
- `dir`: 2D cut plane orientation from 3D outputs ["x","y","z"].
- `sequence`: sequence of plane from - to + in that direction.
- `multifigure`: 1 for multifigure display, 0 for subplots.
- `verbose`: display additional information.
- `density`: density for streamlines.
- `stride`: quiver strides in number of cells.
Right now this can only deal with 2D plots or 3D cuts. Full 3D plots may be supported in the
future.
"""
function plotdata(bd::BATLData, func::AbstractString; dir="x", plotmode="contbar",
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, sequence=1, multifigure=true,
   getrangeonly::Bool=false, levels::Int=0, innermask::Bool=false, verbose::Bool=false,
   density=1.0, stride::Int=10, kwargs...)

   x, w = bd.x, bd.w
   plotmode = split(plotmode)
   vars     = split(func)
   ndim     = bd.head.ndim
   nvar     = length(vars)

   if verbose || getrangeonly
      @info "============ PLOTTING PARAMETERS ==============="
      @info "wnames = $(bd.head.wnames)"
      wmin = Vector{Float64}(undef,nvar)
      wmax = Vector{Float64}(undef,nvar)
      # Display min & max for each variable
      for (ivar,var) in enumerate(vars)
         if occursin(";",var) continue end # skip the vars for streamline
         varIndex_ = findindex(bd, var)
         if ndim == 1
            wmin[ivar] = minimum(w[:,varIndex_])
            wmax[ivar] = maximum(w[:,varIndex_])
         elseif ndim == 2
            wmin[ivar] = minimum(w[:,:,varIndex_])
            wmax[ivar] = maximum(w[:,:,varIndex_])
         end
         println("Min & Max value for $(var) :$(wmin[ivar])",", $(wmax[ivar])")
      end
      if getrangeonly return wmin, wmax end
   end

   ## plot multiple variables with same plotmode
   if length(plotmode) < nvar
      [push!(plotmode, plotmode[i]) for i = 1:nvar-length(plotmode)]
   end

   ## Plot
   if ndim == 1
      for (ivar,var) in enumerate(vars)
         varIndex_ = findindex(bd, var)
         if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end
         if !occursin("scatter", plotmode[ivar])
            plot(x, w[:,varIndex_]; kwargs...)
         else
            scatter(x, w[:,varIndex_]; kwargs...)
         end
         if occursin("grid", plotmode[ivar])
            grid(true)
         end
         xlabel("x"); ylabel("$(var)")
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" bd.head.it bd.head.time
         at = matplotlib.offsetbox.AnchoredText(str,
            loc="lower left", prop=Dict("size"=>8), frameon=true,
            bbox_to_anchor=(0., 1.), bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
      end
   elseif ndim == 2
      for (ivar,var) in enumerate(vars)
         occursin("over", plotmode[ivar]) && (multifigure = false)
         if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end
         if !occursin(";",var)
            varIndex_ = findindex(bd, var)
         end

         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")

            if plotmode[ivar] ∈ ("contbar", "contbarlog")
               c = contourf(bd, var, levels; plotrange, plotinterval, innermask, kwargs...)
            elseif plotmode[ivar] ∈ ("cont", "contlog")
               c = contour(bd, var, levels; plotrange, plotinterval, innermask, kwargs...)
            elseif plotmode[ivar] ∈ ("surfbar", "surfbarlog")
               c = plot_surface(bd, var; plotrange, plotinterval, innermask, kwargs...)
            end

            occursin("bar", plotmode[ivar]) && colorbar(c; fraction=0.04, pad=0.02)
            occursin("log", plotmode[ivar]) && (c.locator = matplotlib.ticker.LogLocator())
            title(bd.head.wnames[varIndex_])

         elseif plotmode[ivar] ∈ ("trimesh","trisurf","tricont","tristream")
            X = vec(x[:,:,1])
            Y = vec(x[:,:,2])
            W = vec(w[:,:,varIndex_])

            # This needs to be modified!!!
            if !all(isinf.(plotrange))
               xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
                  Y .> plotrange[3] .& Y .< plotrange[4]
               X = X[xyIndex]
               Y = Y[xyIndex]
               W = W[xyIndex]
            end

            if plotmode[ivar] == "trimesh"
               triang = matplotlib.tri.Triangulation(X, Y)
               c = ax.triplot(triang, kwargs...)
            elseif plotmode[ivar] == "trisurf"
               c = ax.plot_trisurf(X, Y, W'; kwargs...)
            elseif plotmode[ivar] == "tricont"
               c = ax.tricontourf(X, Y, W'; kwargs...)
               fig.colorbar(c; ax)
            elseif plotmode[ivar] == "tristream"
               throw(ArgumentError("tristream not yet implemented!"))
            end

            title(bd.head.wnames[varIndex_])

         elseif plotmode[ivar] ∈ ("stream","streamover")
            s = streamplot(bd, var; plotrange, plotinterval, color="w", linewidth=1.0,
               density, kwargs...)
         elseif occursin("quiver", plotmode[ivar])
            q = quiver(bd, var; stride, color="w", kwargs...)
         elseif occursin("grid", plotmode[ivar])
            # This does not take subdomain plot into account!
            X, Y = x[:,:,1], x[:,:,2]
            scatter(X, Y, marker=".", alpha=0.6)
            title("Grid illustration")
         else
            throw(ArgumentError("unknown plot mode: $(plotmode[ivar])"))
         end

         xlabel(bd.head.variables[1]); ylabel(bd.head.variables[2])
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" bd.head.it bd.head.time
         at = matplotlib.offsetbox.AnchoredText(str,
            loc="lower left", prop=Dict("size"=>8), frameon=true,
            bbox_to_anchor=(0., 1.), bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
         # recover status
         occursin("over", plotmode[ivar]) && (multifigure = true)
      end

   else # 2D cut from 3D output; now only for Cartesian output
      X = @view x[:,:,:,1]
      Y = @view x[:,:,:,2]
      Z = @view x[:,:,:,3]
      for (ivar,var) in enumerate(vars)
         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar", "contlog",
               "contbarlog")

            varIndex_ = findindex(bd, var)

            if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end

            W = w[:,:,:,varIndex_]

            if dir == "x"
               cut1 = @view X[sequence,:,:]
               cut2 = @view Y[sequence,:,:]
               W    = @view W[sequence,:,:]
            elseif dir == "y"
               cut1 = @view X[:,sequence,:]
               cut2 = @view Z[:,sequence,:]
               W    = @view W[:,sequence,:]
            elseif dir == "z"
               cut1 = @view X[:,:,sequence]
               cut2 = @view Y[:,:,sequence]
               W    = @view W[:,:,sequence]
            end
         elseif plotmode[ivar] ∈ ("stream","streamover")
            varstream  = split(var,";")
            var1_ = findindex(bd, varstream[1])
            var2_ = findindex(bd, varstream[2])

            v1 = @view w[:,:,:,var1_]
            v2 = @view w[:,:,:,var2_]

            if dir == "x"
               cut1 = @view Y[sequence,:,:]
               cut2 = @view Z[sequence,:,:]
               v1   = v1[sequence,:,:]'
               v2   = v2[sequence,:,:]'
            elseif dir == "y"
               cut1 = @view X[:,sequence,:]
               cut2 = @view Z[:,sequence,:]
               v1   = v1[:,sequence,:]'
               v2   = v2[:,sequence,:]'
            elseif dir == "z"
               cut1 = @view X[:,:,sequence]
               cut2 = @view Y[:,:,sequence]
               v1   = v1[:,:,sequence]'
               v2   = v2[:,:,sequence]'
            end
            cut1, cut2 = cut1', cut2'
         end

         if !all(isinf.(plotrange))
            cut1, cut2, v1, v2 = subsurface(cut1, cut2, v1, v2, plotrange)
         end

         if plotmode[ivar] ∈ ("surf", "surfbar", "surfbarlog", "cont", "contbar", "contlog",
               "contbarlog")
            if levels == 0
               c = ax.contourf(cut1, cut2, W'; kwargs...)
            else
               c = ax.contourf(cut1, cut2, W', levels; kwargs...)
            end
            fig.colorbar(c, ax=ax)
            title(bd.head.wnames[varIndex_])

         elseif plotmode[ivar] ∈ ("stream", "streamover")
            xi = range(cut1[1,1], stop=cut1[1,end],
               step=(cut1[1,end] - cut1[1,1]) / (size(cut1,2) - 1))
            yi = range(cut2[1,1], stop=cut2[end,1],
               step=(cut2[end,1] - cut2[1,1]) / (size(cut2,1) - 1))

            s = streamplot(xi, yi, v1, v2; color="w", linewidth=1.0, density, kwargs...)
         end

         if dir == "x"
            xlabel("y"); ylabel("z")
         elseif dir == "y"
            xlabel("x"); ylabel("z")
         elseif dir == "z"
            xlabel("x"); ylabel("y")
         end

         ax = gca()
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" bd.head.it bd.head.time
         at = matplotlib.offsetbox.AnchoredText(str,
            loc="lower left", prop=Dict("size"=>8), frameon=true,
            bbox_to_anchor=(0., 1.), bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
      end
   end

   return
end


"""
    cutplot(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1,
       levels=20)

2D plane cut contourf of 3D box data.
"""
function cutplot(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1, levels::Int=20)

   x, w = bd.x, bd.w
   varIndex_ = findindex(bd, var)

   X = @view x[:,:,:,1]
   Y = @view x[:,:,:,2]
   Z = @view x[:,:,:,3]

   W = @view w[:,:,:,varIndex_]

   if dir == "x"
      cut1 = @view X[sequence,:,:]
      cut2 = @view Y[sequence,:,:]
      W    = @view W[sequence,:,:]
   elseif dir == "y"
      cut1 = @view X[:,sequence,:]
      cut2 = @view Z[:,sequence,:]
      W    = @view W[:,sequence,:]
   elseif dir == "z"
      cut1 = @view X[:,:,sequence]
      cut2 = @view Y[:,:,sequence]
      W    = @view W[:,:,sequence]
   end

   if !all(isinf.(plotrange))
      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
   end
   if isnothing(ax) ax = plt.gca() end
   c = ax.contourf(cut1, cut2, W, levels)

   title(bd.head.wnames[varIndex_])

   if dir == "x"
      xlabel("y"); ylabel("z")
   elseif dir == "y"
      xlabel("x"); ylabel("z")
   elseif dir == "z"
      xlabel("x"); ylabel("y")
   end

   c
end


"""
    streamslice(data::BATLData, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], dir="x",
       sequence=1; kwargs...)

Plot streamlines on 2D slices of 3D box data. Variable names in `var` string must be
separated with `;`.
"""
function streamslice(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1, kwargs...)

   x,w = bd.x, bd.w
   varstream  = split(var, ";")
   var1_ = findindex(bd, varstream[1])
   var2_ = findindex(bd, varstream[2])

   X = @view x[:,:,:,1]
   Y = @view x[:,:,:,2]
   Z = @view x[:,:,:,3]

   v1 = @view w[:,:,:,var1_]
   v2 = @view w[:,:,:,var2_]

   if dir == "x"
      cut1 = @view X[sequence,:,:]
      cut2 = @view Y[sequence,:,:]
      v1   = v1[sequence,:,:]
      v2   = v2[sequence,:,:]
   elseif dir == "y"
      cut1 = @view X[:,sequence,:]
      cut2 = @view Z[:,sequence,:]
      v1   = v1[:,sequence,:]
      v2   = v2[:,sequence,:]
   elseif dir == "z"
      cut1 = @view X[:,:,sequence]
      cut2 = @view Y[:,:,sequence]
      v1   = v1[:,:,sequence]
      v2   = v2[:,:,sequence]
   end

   if !all(isinf.(plotrange))
      cut1, cut2, v1, v2 = subsurface(cut1, cut2, v1, v2, plotrange)
   end

   xi = range(cut1[1,1], stop=cut1[end,1],
      step = (cut1[end,1] - cut1[1,1]) / (size(cut1,1) - 1))
   yi = range(cut2[1,1], stop=cut2[1,end],
      step = (cut2[1,end] - cut2[1,1]) / (size(cut2,2) - 1))

   if isnothing(ax) ax = plt.gca() end
   s = ax.streamplot(xi, yi, v1', v2'; kwargs...)

   if dir == "x"
      xlabel("y"); ylabel("z")
   elseif dir == "y"
      xlabel("x"); ylabel("z")
   elseif dir == "z"
      xlabel("x"); ylabel("y")
   end

   s
end


"""
    plot(data, var, ax=nothing; kwargs...)

Wrapper over `plot` in matplotlib.
"""
function PyPlot.plot(bd::BATLData, var::AbstractString, ax=nothing; kwargs...)
   x, w = bd.x, bd.w
   varIndex_ = findindex(bd, var)
   if isnothing(ax) ax = plt.gca() end

   c = ax.plot(x, w[:,varIndex_]; kwargs...)
end

"""
    scatter(data, var, ax=nothing; kwargs...)

Wrapper over `scatter` in matplotlib.
"""
function PyPlot.scatter(bd::BATLData, var::AbstractString, ax=nothing; kwargs...)
   x, w = bd.x, bd.w
   varIndex_ = findindex(bd, var)
   if isnothing(ax) ax = plt.gca() end

   c = ax.scatter(x, w[:,varIndex_]; kwargs...)
end

"""
    contour(data, var, levels=0; ax=nothing, plotrange=[-Inf,Inf,-Inf,Inf],
       plotinterval=0.1, innermask=false, kwargs...)

Wrapper over `contour` in matplotlib.
"""
function PyPlot.contour(bd::BATLData, var::AbstractString, levels::Int=0; ax=nothing,
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, innermask=false, kwargs...)

   Xi, Yi, Wi = getdata2d(bd, var, plotrange, plotinterval; innermask)

   if isnothing(ax) ax = plt.gca() end

   if levels != 0
      c = ax.contour(Xi, Yi, Wi, levels; kwargs...)
   else
      c = ax.contour(Xi, Yi, Wi; kwargs...)
   end
end

"""
    contourf(data, var, levels=0; ax=nothing, plotrange=[-Inf,Inf,-Inf,Inf],
       plotinterval=0.1, innermask=false, kwargs...)

Wrapper over `contourf` in matplotlib.
"""
function PyPlot.contourf(bd::BATLData, var::AbstractString, levels::Int=0; ax=nothing,
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, innermask=false, kwargs...)

   Xi, Yi, Wi = getdata2d(bd, var, plotrange, plotinterval; innermask)

   if isnothing(ax) ax = plt.gca() end

   if levels != 0
      c = ax.contourf(Xi, Yi, Wi, levels; kwargs...)
   else
      c = ax.contourf(Xi, Yi, Wi; kwargs...)
   end
end

"""
    tricontourf(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
       kwargs...)

Wrapper over `tricontourf` in matplotlib.
"""
function PyPlot.tricontourf(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs...)

   x, w = bd.x, bd.w
   varIndex_ = findindex(bd, var)

   X = vec(x[:,:,1])
   Y = vec(x[:,:,2])
   W = vec(w[:,:,varIndex_])

   # This needs to be modified!!!
   if !all(isinf.(plotrange))
      xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
         Y .> plotrange[3] .& Y .< plotrange[4]
      X = X[xyIndex]
      Y = Y[xyIndex]
      W = W[xyIndex]
   end
   if isnothing(ax) ax = plt.gca() end

   ax.tricontourf(X, Y, W; kwargs...)
end

"""
    plot_trisurf(data::BATLData, var::String, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf],
       kwargs...)

Wrapper over `plot_trisurf` in matplotlib.
"""
function PyPlot.plot_trisurf(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], kwargs...)

   x, w = bd.x, bd.w
   varIndex_ = findindex(bd, var)

   X = vec(x[:,:,1])
   Y = vec(x[:,:,2])
   W = vec(w[:,:,varIndex_])

   # This needs to be modified!!!
   if !all(isinf.(plotrange))
      xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
         Y .> plotrange[3] .& Y .< plotrange[4]
      X = X[xyIndex]
      Y = Y[xyIndex]
      W = W[xyIndex]
   end
   if isnothing(ax) ax = plt.gca() end

   ax.plot_trisurf(X, Y, W)
end


"""
    plot_surface(data, var; plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
       innermask=false, kwargs...)

Wrapper over `plot_surface` in matplotlib.
"""
function PyPlot.plot_surface(bd::BATLData, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, innermask=false, kwargs...)

   xi, yi, Wi = getdata2d(bd, var, plotrange, plotinterval; innermask)
   Xi, Yi = meshgrid(xi, yi)

   plot_surface(Xi, Yi, Wi; kwargs...)
end


"""
    pcolormesh(data, var, levels=0; ax=nothing, plotrange=[-Inf,Inf,-Inf,Inf],
       plotinterval=0.1, innermask=false, kwargs...)

Wrapper over `pcolormesh` in matplotlib.
"""
function PyPlot.pcolormesh(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, innermask=false, kwargs...)

   xi, yi, Wi = getdata2d(bd, var, plotrange, plotinterval; innermask)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(xi, yi, Wi; kwargs...)
end


"""
    tripcolor(data, var, levels=0; ax=nothing, plotrange=[-Inf,Inf,-Inf,Inf],
       plotinterval=0.1, innermask=false, kwargs...)

Wrapper over `tripcolor` in matplotlib.
"""
function PyPlot.tripcolor(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], innermask=false, kwargs...)

   x, w, ndim = bd.x, bd.w, bd.head.ndim

   varIndex_ = findindex(bd, var)

   X = @view x[:,:,1]
   Y = @view x[:,:,2]
   W = vec(w[:,:,varIndex_])

   if any(abs.(plotrange) .== Inf)
      if plotrange[1] == -Inf plotrange[1] = minimum(X) end
      if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
      if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
      if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
   end

   triang = matplotlib.tri.Triangulation(vec(X), vec(Y))

   # Mask off unwanted triangles at the inner boundary.
   if innermask
      varIndex_ = findlast(x->x=="rbody", bd.head.variables)
      isnothing(varIndex_) && error("rbody not found in file header parameters!")
      ParamIndex_ = varIndex_ - ndim - bd.head.nw
      r2 = bd.head.eqpar[ParamIndex_]^2

      ids = triang.triangles .+ Int32(1)
      mask = Vector{Bool}(undef, size(ids, 1))
      for i in axes(ids, 1)
         xmean = sum(@view X[ids[i,:]]) / 3
         ymean = sum(@view Y[ids[i,:]]) / 3
         xmean^2 + ymean^2 < r2 && (mask[i] = true)
      end

      triang.set_mask(mask)
   end


   if isnothing(ax) ax = plt.gca() end

   c = ax.tripcolor(triang, W; kwargs...)

   ax.set_xlim(plotrange[1], plotrange[2])
   ax.set_ylim(plotrange[3], plotrange[4])

   c
end


"""
    streamplot(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
       kwargs...)

Wrapper over `streamplot` in matplotlib .
"""
function PyPlot.streamplot(bd::BATLData, var::AbstractString, ax=nothing;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs...)

   x, w = bd.x, bd.w
   varstream  = split(var,";")
   wnames = lowercase.(bd.head.wnames)
   var1_ = findfirst(x->x==lowercase(varstream[1]), wnames)
   var2_ = findfirst(x->x==lowercase(varstream[2]), wnames)

   if bd.head.gencoord # generalized coordinates
      X, Y = vec(x[:,:,1]), vec(x[:,:,2])
      if any(isinf.(plotrange))
         if plotrange[1] == -Inf plotrange[1] = minimum(X) end
	      if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
         if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
	      if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
      end

      # Create grid values first.
      xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
      yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

      # Is there a triangulation method in Julia?
      tr = matplotlib.tri.Triangulation(X, Y)
      Xi, Yi = meshgrid(xi, yi)

      interpolator = matplotlib.tri.LinearTriInterpolator(tr, w[:,1,var1_])
      v1 = interpolator(Xi, Yi)

      interpolator = matplotlib.tri.LinearTriInterpolator(tr, w[:,1,var2_])
      v2 = interpolator(Xi, Yi)

   else # Cartesian coordinates
      xrange = range(x[1,1,1], x[end,1,1], size(x,1))
      yrange = range(x[1,1,2], x[1,end,2], size(x,2))
      if all(isinf.(plotrange))
         xi = xrange
         yi = yrange
         v1 = w[:,:,var1_]'
         v2 = w[:,:,var2_]'
      else
	      if plotrange[1] == -Inf plotrange[1] = xrange[1] end
	      if plotrange[2] ==  Inf plotrange[2] = xrange[end] end
         if plotrange[3] == -Inf plotrange[3] = yrange[1] end
         if plotrange[4] ==  Inf plotrange[4] = yrange[end] end

         w1, w2 = @views w[:,:,var1_], w[:,:,var2_]
         xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
         yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
         interp1 = cubic_spline_interpolation((xrange, yrange), w1)
         v1 = [interp1(i, j) for j in yi, i in xi]
   
         interp2 = cubic_spline_interpolation((xrange, yrange), w2)
         v2 = [interp2(i, j) for j in yi, i in xi]
      end
   end
   if isnothing(ax) ax = plt.gca() end

   ax.streamplot(xi, yi, v1, v2; kwargs...)
end

"""
    quiver(data, var, ax=nothing; stride=10, kwargs...)

Wrapper over `quiver` in matplotlib. Only supports Cartesian grid for now.
"""
function PyPlot.quiver(bd::BATLData, var::AbstractString, ax=nothing;
   stride::Integer=10, kwargs...)
   x, w = bd.x, bd.w
   VarQuiver  = split(var, ";")
   var1_ = findindex(bd, VarQuiver[1])
   var2_ = findindex(bd, VarQuiver[2])

   @views X, Y = x[:,1,1], x[1,:,2]
   @views v1, v2 = w[:,:,var1_]', w[:,:,var2_]'

   @views Xq, Yq = X[1:stride:end], Y[1:stride:end]
   v1q = @view v1[1:stride:end, 1:stride:end]
   v2q = @view v2[1:stride:end, 1:stride:end]

   if isnothing(ax) ax = plt.gca() end

   ax.quiver(Xq, Yq, v1q, v2q; kwargs...)
end

"Set colorbar norm and ticks."
function set_colorbar(colorscale, vmin, vmax, data=[1.0])
   if colorscale == :linear || any(<(0), data)
      colorscale == :log && @warn "Nonpositive data detected: use linear scale instead!"
      v1 = isinf(vmin) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
      nticks = 7
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(v1, v2)
      cnorm = matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
      ticks = range(v1, v2, length=nticks)
   else # logarithmic
      datapositive = data[data .> 0.0]
      v1 = isinf(vmin) ? minimum(datapositive) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax

      cnorm = matplotlib.colors.LogNorm(vmin=v1, vmax=v2)
      ticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
   end

   cnorm, ticks
end