# Plotting functionalities.

using PyPlot
using Dierckx: Spline2D

export plotdata, plotlogdata, animatedata, plot, scatter, contour, contourf,
	   plot_surface, tricontourf, plot_trisurf, streamplot, streamslice, cutplot

import PyPlot: plot, scatter, contour, contourf, plot_surface, tricontourf,
	   plot_trisurf, streamplot

include("utility.jl")

"""
    plotlogdata(data, head, func; plotmode="line")

Plot information from log file.
# Input arguments
- `data::Array`: output data.
- `head::NamedTuple`: header info.
- `func::String`: variables for plotting.
- `plotmode::String`: type of plotting ["line","scatter"].
"""
function plotlogdata(data, head::NamedTuple, func::AbstractString;
   plotmode="line")

   vars     = split(func)
   plotmode = split(plotmode)

   for (ivar, var) in enumerate(vars)
      VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(head.variables))
      isnothing(VarIndex_) && error("$(var) not found in file header variables!")

      figure()
      if plotmode[ivar] == "line"
         plot(data[1,:],data[VarIndex_,:])
      elseif plotmode[ivar] == "scatter"
         scatter(data[1,:],data[VarIndex_,:])
      else
         error("unknown plot mode for plotlogdata!")
      end
      xlabel(head.variables[1])
      ylabel(head.variables[VarIndex_])
      title("log file data")
   end

end


"""
    plotdata(data, func, (...))

Plot the variable from SWMF output.

`plotdata(data, "p", plotmode="contbar")`

`plotdata(data, "p", plotmode="grid")`

`plotdata(data, func, plotmode="trimesh",plotrange=plotrange, plotinterval=0.2)`

# Input arguments
- `data::Data`: output data.
- `func::String`: variables for plotting.
- `plotmode::String`: (optional) type of plotting ["cont","contbar"]...
- `plotrange::Vector`: (optional) range of plotting.
- `plotinterval`: (optional) interval for interpolation.
- `level`: (optional) level of contour.
- `cut`: (optional) select 2D cut plane from 3D outputs ["x","y","z"].
- `cutPlaneIndex`: (optional)
- `multifigure`: (optional) 1 for multifigure display, 0 for subplots.
- `verbose`: (optional) display additional information.
- `density`: (optional) density for streamlines.
Right now this can only deal with 2D plots or 3D cuts. Full 3D plots may be
supported in the future.
"""
function plotdata(data::Data, func::AbstractString; cut="", plotmode="contbar",
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, cutPlaneIndex=1,
   multifigure=true, getrangeOnly=false, level=0, verbose=false, density=1.0)

   x, w = data.x, data.w
   plotmode = split(plotmode)
   vars     = split(func)
   ndim     = data.head.ndim
   nvar     = length(vars)

   if verbose || getrangeOnly
      @info "============ PLOTTING PARAMETERS ==============="
      @info "wnames = $(data.head.wnames)"
      wmin = Vector{Float64}(undef,nvar)
      wmax = Vector{Float64}(undef,nvar)
      # Display min & max for each variable
      for (ivar,var) in enumerate(vars)
         if occursin(";",var) continue end # skip the vars for streamline
         VarIndex_ = findindex(data, var)
         if ndim == 1
            wmin[ivar] = minimum(w[:,VarIndex_])
            wmax[ivar] = maximum(w[:,VarIndex_])
         elseif ndim == 2
            wmin[ivar] = minimum(w[:,:,VarIndex_])
            wmax[ivar] = maximum(w[:,:,VarIndex_])
         end
         println("Min & Max value for $(var) :$(wmin[ivar])",", $(wmax[ivar])")
      end
      if getrangeOnly return wmin, wmax end
   end

   ## plot multiple variables with same plotmode
   if length(plotmode) < nvar
      [push!(plotmode, plotmode[i]) for i = 1:nvar-length(plotmode)]
   end

   ## Plot
   if ndim == 1
      for (ivar,var) in enumerate(vars)
         VarIndex_ = findindex(data, var)
         if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end
         if !occursin("scatter",plotmode[ivar])
            plot(x,w[:,VarIndex_])
         else
            scatter(x,w[:,VarIndex_])
         end
         if occursin("grid",plotmode[ivar])
            grid(true)
         end
         xlabel("x"); ylabel("$(var)")
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" data.head.it data.head.time
         at = matplotlib.offsetbox.AnchoredText(str,
         loc="lower left", prop=Dict("size"=>8), frameon=true,
         bbox_to_anchor=(0., 1.),
         bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
      end
   elseif ndim == 2
      for (ivar,var) in enumerate(vars)
         occursin("over", plotmode[ivar]) && (multifigure = false)
         if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end
         if !occursin(";",var)
            VarIndex_ = findindex(data, var)
         end

         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")

            Xi, Yi, Wi = getdata(data, var, plotrange, plotinterval, 2)

            # More robust method needed!
            if plotmode[ivar] ∈ ["contbar", "contbarlog"]
               if level == 0
                  c = contourf(Xi, Yi, Wi)
               else
                  c = contourf(Xi, Yi, Wi, level)
               end
            elseif plotmode[ivar] ∈ ["cont", "contlog"]
               c = contour(Xi, Yi, Wi)
            elseif plotmode[ivar] ∈ ["surfbar", "surfbarlog"]
               c = plot_surface(Xi, Yi, Wi)
            end

            occursin("bar", plotmode[ivar]) && colorbar()
            occursin("log", plotmode[ivar]) &&
            ( c.locator = matplotlib.ticker.LogLocator() )
            title(data.head.wnames[VarIndex_])

         elseif plotmode[ivar] ∈ ("trimesh","trisurf","tricont","tristream")
            X = vec(x[:,:,1])
            Y = vec(x[:,:,2])
            W = vec(w[:,:,VarIndex_])

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
               c = ax.triplot(triang)
            elseif plotmode[ivar] == "trisurf"
               c = ax.plot_trisurf(X, Y, W)
            elseif plotmode[ivar] == "tricont"
               c = ax.tricontourf(X, Y, W)
               fig.colorbar(c,ax=ax)
            elseif plotmode[ivar] == "tristream"
               error("not yet implemented!")
            end

            title(data.head.wnames[VarIndex_])

         elseif plotmode[ivar] ∈ ("stream","streamover")
            VarStream  = split(var,";")
            VarIndex1_ = findindex(data, VarStream[1])
            VarIndex2_ = findindex(data, VarStream[2])

            if data.head.gencoord # Generalized coordinates
               xrange = @view x[:,:,1]
               yrange = @view x[:,:,2]
               if any(isinf.(plotrange))
                  if plotrange[1] == -Inf plotrange[1] = minimum(xrange) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(xrange) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(yrange) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(yrange) end
               end

               # Create grid values first.
               xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
               yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

               # The PyCall here can be potentially replaced with Spline2D.
               # Perform linear interpolation of the data (x,y) on grid(xi,yi)
               triang = matplotlib.tri.Triangulation(xrange[:], yrange[:])

               Xi, Yi = meshgrid(xi, yi)

               W = w[:,1,VarIndex1_]
               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               v1 = interpolator(Xi, Yi)

               W = w[:,1,VarIndex2_]
               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               v2 = interpolator(Xi, Yi)

            else # Cartesian coordinates
               xrange = @view x[:,1,1]
               yrange = @view x[1,:,2]
               if all(isinf.(plotrange))
                  Xi, Yi = meshgrid(xrange, yrange)
                  v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'
               else
                  if plotrange[1] == -Inf plotrange[1] = minimum(xrange) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(xrange) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(yrange) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(yrange) end

                  w1, w2 = w[:,:,VarIndex1_], w[:,:,VarIndex2_]

                  xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
                  yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

                  Xi, Yi = meshgrid(xi, yi)

                  spline = Spline2D(xrange, yrange, w1)
                  v1 = spline(Xi[:], Yi[:])
                  v1 = reshape(v1, size(Xi))

                  spline = Spline2D(xrange, yrange, w2)
                  v2 = spline(Xi[:], Yi[:])
                  v2 = reshape(v2, size(Xi))
               end
            end

            s = streamplot(Xi,Yi,v1,v2; color="w", linewidth=1.0, density)

         elseif occursin("quiver", plotmode[ivar])
            VarQuiver  = split(var,";")
            VarIndex1_ = findindex(data, VarQuiver[1])
            VarIndex2_ = findindex(data, VarQuiver[2])

            X, Y = x[:,1,1], x[1,:,2]
            v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'

            q = quiver(X,Y,v1,v2,color="w")

         elseif occursin("grid", plotmode[ivar])
            # This does not take subdomain plot into account!
            X, Y = x[:,:,1], x[:,:,2]
            scatter(X,Y,marker=".",alpha=0.6)
            title("Grid illustration")
         else
            error("unknown plot mode: $(plotmode[ivar])")
         end

         xlabel(data.head.variables[1]); ylabel(data.head.variables[2])
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" data.head.it data.head.time
         at = matplotlib.offsetbox.AnchoredText(str,
         loc="lower left", prop=Dict("size"=>8), frameon=true,
         bbox_to_anchor=(0., 1.),
         bbox_transform=ax.transAxes)
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
         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")

            VarIndex_ = findindex(data, var)

            if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end

            W = w[:,:,:,VarIndex_]

            if cut ∈ ("x","")
               cut1 = @view X[cutPlaneIndex,:,:]
               cut2 = @view Y[cutPlaneIndex,:,:]
               W    = @view W[cutPlaneIndex,:,:]
            elseif cut ==  "y"
               cut1 = @view X[:,cutPlaneIndex,:]
               cut2 = @view Z[:,cutPlaneIndex,:]
               W    = @view W[:,cutPlaneIndex,:]
            elseif cut == "z"
               cut1 = @view X[:,:,cutPlaneIndex]
               cut2 = @view Y[:,:,cutPlaneIndex]
               W    = @view W[:,:,cutPlaneIndex]
            end
         elseif plotmode[ivar] ∈ ("stream","streamover")
            varStream  = split(var,";")
            VarIndex1_ = findindex(data, varStream[1])
            VarIndex2_ = findindex(data, varStream[2])

            v1 = @view w[:,:,:,VarIndex1_]
            v2 = @view w[:,:,:,VarIndex2_]

            if cut ∈ ("x","")
               cut1 = @view Y[cutPlaneIndex,:,:]
               cut2 = @view Z[cutPlaneIndex,:,:]
               v1   = v1[cutPlaneIndex,:,:]'
               v2   = v2[cutPlaneIndex,:,:]'
            elseif cut ==  "y"
               cut1 = @view X[:,cutPlaneIndex,:]
               cut2 = @view Z[:,cutPlaneIndex,:]
               v1   = v1[:,cutPlaneIndex,:]'
               v2   = v2[:,cutPlaneIndex,:]'
            elseif cut == "z"
               cut1 = @view X[:,:,cutPlaneIndex]
               cut2 = @view Y[:,:,cutPlaneIndex]
               v1   = v1[:,:,cutPlaneIndex]'
               v2   = v2[:,:,cutPlaneIndex]'
            end
            cut1, cut2 = cut1', cut2'
         end

         if !all(isinf.(plotrange))
            cut1, cut2, v1, v2 = subsurface(cut1, cut2, v1, v2, plotrange)
         end

         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")
            if level == 0
               c = ax.contourf(cut1, cut2, W)
            else
               c = ax.contourf(cut1, cut2, W, level)
            end
            fig.colorbar(c, ax=ax)
            title(data.head.wnames[VarIndex_])

         elseif plotmode[ivar] ∈ ("stream","streamover")
            # Surprisingly, some box outputs do not have equal spaces???
            #xi = range(cut1[1,1], stop=cut1[1,end], length=size(cut1)[2])
            #yi = range(cut2[1,1], stop=cut2[end,1], length=size(cut2)[1])

            xi = range(cut1[1,1], stop=cut1[1,end],
            step=(cut1[1,end] - cut1[1,1]) / (size(cut1,2) - 1))
            yi = range(cut2[1,1], stop=cut2[end,1],
            step=(cut2[end,1] - cut2[1,1]) / (size(cut2,1) - 1))

            Xi, Yi = meshgrid(xi, yi)

            s = streamplot(Xi,Yi,v1,v2; color="w", linewidth=1.0, density)
         end

         if cut == "x"
            xlabel("y"); ylabel("z")
         elseif cut == "y"
            xlabel("x"); ylabel("z")
         elseif cut == "z"
            xlabel("x"); ylabel("y")
         end

         ax = gca()
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" data.head.it data.head.time
         at = matplotlib.offsetbox.AnchoredText(str,
         loc="lower left", prop=Dict("size"=>8), frameon=true,
         bbox_to_anchor=(0., 1.),
         bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
      end
   end

end


"""
    cutplot(data, var; plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',
       plotinterval=0.1, density=1.0, cutPlaneIndex=1,level=20)

2D plane cut contourf of 3D box data.
"""
function cutplot(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',
   plotinterval=0.1, density=1.0, cutPlaneIndex=1,level=20)

   x, w = data.x, data.w
   VarIndex_ = findindex(data, var)

   X = @view x[:,:,:,1]
   Y = @view x[:,:,:,2]
   Z = @view x[:,:,:,3]

   W = w[:,:,:,VarIndex_]

   if cut ∈ ('x',' ')
      cut1 = @view X[cutPlaneIndex,:,:]
      cut2 = @view Y[cutPlaneIndex,:,:]
      W    = @view W[cutPlaneIndex,:,:]
   elseif cut ==  'y'
      cut1 = @view X[:,cutPlaneIndex,:]
      cut2 = @view Z[:,cutPlaneIndex,:]
      W    = @view W[:,cutPlaneIndex,:]
   elseif cut == 'z'
      cut1 = @view X[:,:,cutPlaneIndex]
      cut2 = @view Y[:,:,cutPlaneIndex]
      W    = @view W[:,:,cutPlaneIndex]
   end

   if !all(isinf.(plotrange))
      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
   end

   c = contourf(cut1, cut2, W, level)

   title(data.head.wnames[VarIndex_])

   if cut == "x"
      xlabel("y"); ylabel("z")
   elseif cut == "y"
      xlabel("x"); ylabel("z")
   elseif cut == "z"
      xlabel("x"); ylabel("y")
   end

   return c::PyCall.PyObject
end


"""
    streamslice(data::Data, var::String;
       plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',
       plotinterval=0.1, density=1.0, cutPlaneIndex=1, color="w", linewidth=1.0)

Plot streamlines on 2D slices of 3D box data. Variable string must be separated
with `;`. Tranposes aree needed because of `meshgrid` and `ndgrid` conversion.
"""
function streamslice(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], cut=' ', cutPlaneIndex=1,
   plotinterval=0.1, kwargs...)

   x,w = data.x, data.w
   varStream  = split(var,";")
   VarIndex1_ = findindex(data, varStream[1])
   VarIndex2_ = findindex(data, varStream[2])

   X = @view x[:,:,:,1]
   Y = @view x[:,:,:,2]
   Z = @view x[:,:,:,3]

   v1 = @view w[:,:,:,VarIndex1_]
   v2 = @view w[:,:,:,VarIndex2_]

   if cut ∈ ('x',' ')
      cut1 = @view X[cutPlaneIndex,:,:]
      cut2 = @view Y[cutPlaneIndex,:,:]
      v1   = v1[cutPlaneIndex,:,:]
      v2   = v2[cutPlaneIndex,:,:]
   elseif cut ==  'y'
      cut1 = @view X[:,cutPlaneIndex,:]
      cut2 = @view Z[:,cutPlaneIndex,:]
      v1   = v1[:,cutPlaneIndex,:]
      v2   = v2[:,cutPlaneIndex,:]
   elseif cut == 'z'
      cut1 = @view X[:,:,cutPlaneIndex]
      cut2 = @view Y[:,:,cutPlaneIndex]
      v1   = v1[:,:,cutPlaneIndex]
      v2   = v2[:,:,cutPlaneIndex]
   end

   if !all(isinf.(plotrange))
      cut1, cut2, v1, v2 = subsurface(cut1, cut2, v1, v2, plotrange)
   end

   xi = range(cut1[1,1], stop=cut1[end,1],
   step = (cut1[end,1] - cut1[1,1]) / (size(cut1,1) - 1))
   yi = range(cut2[1,1], stop=cut2[1,end],
   step = (cut2[1,end] - cut2[1,1]) / (size(cut2,2) - 1))

   Xi, Yi = meshgrid(xi, yi)

   s = streamplot(Xi,Yi,v1',v2'; kwargs...)

   if cut == 'x'
      xlabel("y"); ylabel("z")
   elseif cut == 'y'
      xlabel("x"); ylabel("z")
   elseif cut == 'z'
      xlabel("x"); ylabel("y")
   end
   s
end


"""
    plot(data, var; kwargs...)

Wrapper over the plot function in matplotlib.
"""
function plot(data::Data, var::AbstractString; kwargs...)
   x, w = data.x, data.w
   VarIndex_ = findindex(data, var)

   c = plot(x, w[:,VarIndex_]; kwargs...)

   return c::Vector{PyCall.PyObject}
end

"""
    scatter(data, var; kwargs...)

Wrapper over the scatter function in matplotlib.
"""
function scatter(data::Data, var::AbstractString; kwargs...)
   x, w = data.x, data.w
   VarIndex_ = findindex(data, var)

   c = scatter(x, w[:,VarIndex_]; kwargs...)

   return c::Vector{PyCall.PyObject}
end

"""
    contour(data, var, levels=0; plotrange, plotinterval, kwargs...)

Wrapper over the contour function in matplotlib.
"""
function contour(data::Data, var::AbstractString, levels=0;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs...)

   Xi, Yi, Wi = getdata(data, var, plotrange, plotinterval)

   if levels != 0
      c = plt.contour(Xi, Yi, Wi, levels; kwargs...)
   else
      c = plt.contour(Xi, Yi, Wi; kwargs...)
   end
end

"""
    contourf(data, var, levels=0; plotrange, plotinterval, kwargs...)

Wrapper over the contourf function in matplotlib.
"""
function contourf(data::Data, var::AbstractString, levels=0;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs...)

   Xi, Yi, Wi = getdata(data, var, plotrange, plotinterval)

   if levels != 0
      c = plt.contourf(Xi, Yi, Wi, levels; kwargs...)
   else
      c = plt.contourf(Xi, Yi, Wi; kwargs...)
   end
end

"""
    tricontourf(data, var; plotrange, plotinterval, kwargs...)

Wrapper over the tricontourf function in matplotlib.
"""
function tricontourf(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs...)

   x, w = data.x, data.w
   VarIndex_ = findindex(data, var)

   X = vec(x[:,:,1])
   Y = vec(x[:,:,2])
   W = vec(w[:,:,VarIndex_])

   # This needs to be modified!!!
   if !all(isinf.(plotrange))
      xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
         Y .> plotrange[3] .& Y .< plotrange[4]
      X = X[xyIndex]
      Y = Y[xyIndex]
      W = W[xyIndex]
   end

   tricontourf(X, Y, W)
end

"""
    plot_trisurf(data::Data, var::String;
       plotrange::Vector{Float64}=[-Inf,Inf,-Inf,Inf], kwargs::Dict=Dict())

Wrapper over the plot_trisurf function in matplotlib.
"""
function plot_trisurf(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], kwargs...)

   x, w = data.x, data.w
   VarIndex_ = findindex(data, var)

   X = vec(x[:,:,1])
   Y = vec(x[:,:,2])
   W = vec(w[:,:,VarIndex_])

   # This needs to be modified!!!
   if !all(isinf.(plotrange))
      xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
         Y .> plotrange[3] .& Y .< plotrange[4]
      X = X[xyIndex]
      Y = Y[xyIndex]
      W = W[xyIndex]
   end

   plot_trisurf(X, Y, W)
end


"""
    plot_surface(data, var; plotrange, plotinterval, kwargs...)

Wrapper over the plot_surface function in matplotlib.
"""
function plot_surface(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs=Dict())

   Xi, Yi, Wi = getdata(data, var, plotrange, plotinterval, 2)

   plot_surface(Xi, Yi, Wi; kwargs...)
end

"""
    streamplot(data, var; plotrange, plotinterval=0.1, kwargs...)

Wrapper over the streamplot function in matplotlib. Streamplot does not have
**kwargs in the API, but it supports `density`, `color`, and some other
keywords.
"""
function streamplot(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs...)

   x, w = data.x, data.w
   VarStream  = split(var,";")
   wnames = lowercase.(data.head.wnames)
   VarIndex1_ = findfirst(x->x==lowercase(VarStream[1]), wnames)
   VarIndex2_ = findfirst(x->x==lowercase(VarStream[2]), wnames)

   if data.head.gencoord # Generalized coordinates
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

      interpolator = matplotlib.tri.LinearTriInterpolator(tr, w[:,1,VarIndex1_])
      v1 = interpolator(Xi, Yi)

      interpolator = matplotlib.tri.LinearTriInterpolator(tr, w[:,1,VarIndex2_])
      v2 = interpolator(Xi, Yi)

   else # Cartesian coordinates
      xrange = @view x[:,1,1]
      yrange = @view x[1,:,2]
      if all(isinf.(plotrange))
         Xi, Yi = x[:,:,1]', x[:,:,2]'
         v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'
      else
	      if plotrange[1] == -Inf plotrange[1] = minimum(xrange) end
	      if plotrange[2] ==  Inf plotrange[2] = maximum(xrange) end
         if plotrange[3] == -Inf plotrange[3] = minimum(yrange) end
         if plotrange[4] ==  Inf plotrange[4] = maximum(yrange) end

	      w1, w2 = w[:,:,VarIndex1_], w[:,:,VarIndex2_]

         xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
         yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

         Xi, Yi = meshgrid(xi, yi)

         spline = Spline2D(xrange, yrange, w1)
	      v1 = spline(Xi[:], Yi[:])
	      v1 = reshape(v1, size(Xi))

	      spline = Spline2D(xrange, yrange, w2)
	      v2 = spline(Xi[:], Yi[:])
	      v2 = reshape(v2, size(Xi))
      end
   end

   streamplot(Xi, Yi, v1, v2; kwargs...)
end