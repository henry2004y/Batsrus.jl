# Plotting functionalities.

using PyPlot

export plotdata, plotlogdata, animatedata, get_vars, cutplot,
       plot, scatter, contour, contourf, plot_surface, tricontourf,
       plot_trisurf, streamplot, subvolume, subsurface

import PyPlot.plot, PyPlot.scatter, PyPlot.contour, PyPlot.contourf,
       PyPlot.plot_surface, PyPlot.tricontourf, PyPlot.plot_trisurf,
       PyPlot.streamplot

"""
	plotlogdata(data, head, vars, (plotmode="line", plotrange=[-Inf,Inf]))

Plot information from log file.
# Input arguments
- `data::Data`: original variable data.
- `head::Dict`: file header information.
- `vars::String`: variables for plotting.
- `plotmode::String`: (optional) type of plotting ["line","scatter"].
- `plotrange::Vector`: (optional) range of plotting.
"""
function plotlogdata(data::Data, head::Dict, func::String;
   plotmode::String="line", plotrange::Vector=[-Inf,Inf] )

   vars     = split(func)
   plotmode = split(plotmode)

   for (ivar, var) in enumerate(vars)
      # find the index for variables
      VarIndex_ = findfirst(x->x==var, head[:variables])

      isnothing(VarIndex_) &&
      error("unknown plotting variable $(func[ivar])!")
      figure()
      if plotmode[ivar] == "line"
         plot(data[:,1],data[:,VarIndex_])
      elseif plotmode[ivar] == "scatter"
         scatter(data[:,1],data[:,VarIndex_])
      else
         error("unknown plot mode for plotlogdata!")
      end
      xlabel(head[:variables][1])
      ylabel(head[:variables][VarIndex_])
      title("log file data")
   end

end


"""
	plotdata(data, head, func, (...))

Plot the variable from SWMF output.

`plotdata(data, head, "p", plotmode="contbar")`

`plotdata(data, head, "p", plotmode="grid")`

`plotdata(data, head, func, plotmode="trimesh",plotrange=plotrange,
   plotinterval=0.2)`

# Input arguments
- `data::Data`: original variable data.
- `head::Dict`: file header information.
- `vars::String`: variables for plotting.
- `plotmode::String`: (optional) type of plotting ["cont","contbar"]...
- `plotrange::Vector`: (optional) range of plotting.
- `plotinterval`: (optional) interval for interpolation.
- `level`: (optional) level of contour.
- `density`: (optional) density for streamlines.
- `cut`: (optional) select 2D cut plane from 3D outputs ["x","y","z"].
- `cutPlaneIndex`: (optional)
- `multifigure`: (optional) 1 for multifigure display, 0 for subplots.
- `verbose`: (optional) display additional information.
Right now this can only deal with 2D plots or 3D cuts. Full 3D plots may be
supported in the future.
"""
function plotdata(data::Data, head::Dict, func::String; cut::String="",
   plotmode="contbar", plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
   density=1.0, cutPlaneIndex=1, multifigure=true, getrangeOnly=false, level=0,
   verbose=false)

   x,w = data.x, data.w
   plotmode = split(plotmode)
   vars     = split(func)
   ndim     = head[:ndim]
   nvar     = length(vars)

   if verbose || getrangeOnly
      println("============ PLOTTING PARAMETERS ===============")
      println("wnames = $(head[:wnames])")
      println("================================================")
      wmin = Vector{Float64}(undef,nvar)
      wmax = Vector{Float64}(undef,nvar)
      # Display min & max for each variable
      for (ivar,var) in enumerate(vars)
         if occursin(";",var) continue end # skip the vars for streamline
         VarIndex_ = findfirst(x->x==lowercase(var),
         lowercase.(head[:wnames]))
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
         VarIndex_ = findfirst(x->x==lowercase(var),
            lowercase.(head[:wnames]))
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
         str = @sprintf "it=%d, time=%4.2f" head[:it] head[:time]
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
            VarIndex_ = findfirst(x->x==lowercase(var),
               lowercase.(head[:wnames]))
            isempty(VarIndex_) &&
            error("$(var) not found in file header variables!")
         end

         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")

            xi, yi, wi = getdata(data, head, var, plotrange, plotinterval)

            # More robust method needed!
            if plotmode[ivar] ∈ ["contbar", "contbarlog"]
			   if level == 0
                  c = contourf(xi, yi, wi)
			   else
                  c = contourf(xi, yi, wi, level)
			   end
            elseif plotmode[ivar] ∈ ["cont", "contlog"]
               c = contour(xi, yi, wi)
            elseif plotmode[ivar] ∈ ["surfbar", "surfbarlog"]
               c = plot_surface(xi, yi, wi)
            end

            occursin("bar", plotmode[ivar]) && colorbar()
            occursin("log", plotmode[ivar]) &&
            ( c.locator = matplotlib.ticker.LogLocator() )
            title(head[:wnames][VarIndex_])

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

            title(head[:wnames][VarIndex_])

         elseif plotmode[ivar] ∈ ("stream","streamover")
            VarStream  = split(var,";")
            VarIndex1_ = findfirst(x->x==lowercase(VarStream[1]),
            lowercase.(head[:wnames]))
            VarIndex2_ = findfirst(x->x==lowercase(VarStream[2]),
            lowercase.(head[:wnames]))

            if head[:gencoord] # Generalized coordinates
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

               # The PyCall here can be potentially replaced with Spline2D.
               # Perform linear interpolation of the data (x,y) on grid(xi,yi)
               triang = matplotlib.tri.Triangulation(X,Y)
               Xi = [y for x in xi, y in yi]
          	   Yi = [x for x in xi, y in yi]

               W = w[:,1,VarIndex1_]
               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               v1 = interpolator(Xi, Yi)

               W = w[:,1,VarIndex2_]
               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               v2 = interpolator(Xi, Yi)

            else # Cartesian coordinates
               X = x[:,1,1]
               Y = x[1,:,2]
               if all(isinf.(plotrange))
                  Xi, Yi = X, Y
                  v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'
               else
                  if plotrange[1] == -Inf plotrange[1] = minimum(X) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

                  w1, w2 = w[:,:,VarIndex1_], w[:,:,VarIndex2_]

                  xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
                  yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

                  Xi = [i for i in xi, j in yi]
                  Yi = [j for i in xi, j in yi]

                  spline = Spline2D(X, Y, w1)
                  v1 = spline(Xi[:], Yi[:])
                  v1 = reshape(v1, size(Xi))'

                  spline = Spline2D(X, Y, w2)
                  v2 = spline(Xi[:], Yi[:])
                  v2 = reshape(v2, size(Xi))'
               end
            end

            s = streamplot(Xi,Yi,v1,v2,color="w",linewidth=1.0,density=density)

         elseif occursin("quiver", plotmode[ivar])
            VarQuiver  = split(var,";")
            VarIndex1_ = findfirst(x->x==lowercase(VarQuiver[1]),
            lowercase.(head[:wnames]))
            VarIndex2_ = findfirst(x->x==lowercase(VarQuiver[2]),
            lowercase.(head[:wnames]))

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

         xlabel(head[:variables][1]); ylabel(head[:variables][2])
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" head[:it] head[:time]
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
            VarIndex_ = findfirst(x->x==lowercase(var),
            lowercase.(head[:wnames]))
            isempty(VarIndex_) && error("$(var) not found in header variables!")

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
            VarStream  = split(var,";")
            VarIndex1_ = findfirst(x->x==lowercase(VarStream[1]),
            lowercase.(head[:wnames]))
            VarIndex2_ = findfirst(x->x==lowercase(VarStream[2]),
            lowercase.(head[:wnames]))
            (isempty(VarIndex1_) || isempty(VarIndex2_)) &&
            error("$(VarStream) not found in header variables!")

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
            cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
         end

         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")
            if level == 0
               c = ax.contourf(cut1, cut2, W)
            else
               c = ax.contourf(cut1, cut2, W, level)
            end
            fig.colorbar(c, ax=ax)
            title(head[:wnames][VarIndex_])

         elseif plotmode[ivar] ∈ ("stream","streamover")
            # Surprisingly, some box outputs do not have equal spaces???
            #xi = range(cut1[1,1], stop=cut1[1,end], length=size(cut1)[2])
            #yi = range(cut2[1,1], stop=cut2[end,1], length=size(cut2)[1])

            xi = range(cut1[1,1], stop=cut1[1,end],
            step=(cut1[1,end]-cut1[1,1])/(size(cut1,2)-1))
            yi = range(cut2[1,1], stop=cut2[end,1],
            step=(cut2[end,1]-cut2[1,1])/(size(cut2,1)-1))

            Xi = [i for j in yi, i in xi]
            Yi = [j for j in yi, i in xi]

            s = streamplot(Xi,Yi,v1,v2,color="w",linewidth=1.0,density=density)
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
         str = @sprintf "it=%d, time=%4.2f" head[:it] head[:time]
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
	subsurface(x, y, data, limits)

Extract subset of 2D surface dataset.
This is a simplified version of subvolume.
"""
function subsurface(x, y, data, limits)

   if length(limits)!=4
      @error "Reduction must be [xmin xmax ymin ymax]"
   end

   if limits[1] > limits[2] || limits[3] > limits[4]
      @error "subsurface:InvalidReductionXRange"
   end

   sz = size(data)

   hx = x[:,1]
   hy = y[1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])

   newdata = subdata(data, xind, yind, sz)

   newx = x[xind, yind]
   newy = y[xind, yind]

   return newx, newy, newdata
end

"""
	subvolume(x, y, z, data, limits)

Extract subset of 3D dataset in ndgrid format.
"""
function subvolume(x, y, z, data, limits)

   if length(limits)!=4
      @error "Reduction must be [xmin xmax ymin ymax]"
   end

   if limits[1] > limits[2] || limits[3] > limits[4] || limits[5] > limits[6]
      @error "subvolume:InvalidReductionXRange"
   end

   sz = size(data)

   hx = x[:,1,1]
   hy = y[1,:,1]
   hz = z[1,1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[5]) limits[5] = maximum(hz) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end
   if isinf(limits[6]) limits[6] = maximum(hz) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])
   zind = findall(limits[5] .≤ hy .≤ limits[6])

   newdata = subdata(data, xind, yind, zind, sz)

   newx = x[xind, yind, zind]
   newy = y[xind, yind, zind]
   newz = z[xind, yind, zind]

   return newx, newy, newz, newdata
end

"""
	subdata(data, xind, yind, sz)
	subdata(data, xind, yind, zind, sz)

Return the sliced data based on indexes.
"""
function subdata(data, xind::Vector{Int}, yind::Vector{Int}, sz::Tuple{Int,Int})
   newdata = data[xind,yind]
   newsz = size(newdata)

   if length(sz) > 2
      newdata = reshape(newdata, (newsz[1:2]..., sz[3:end]))
   end

   return newdata
end

function subdata(data, xind::Vector{Int}, yind::Vector{Int}, zind::Vector{Int},
   sz::Tuple{Int,Int,Int})

   newdata = data[xind,yind,zind]
   newsz = size(newdata)

   if length(sz) > 3
      newdata = reshape(newdata, (newsz[1:3]..., sz[4:end]))
   end

   return newdata
end


"""
	cutplot(data, head, var; plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',
		plotinterval=0.1, density=1.0, cutPlaneIndex=1,level=20)

2D plane cut of 3D box data.
"""
function cutplot(data::Data, head::Dict, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',
   plotinterval=0.1, density=1.0, cutPlaneIndex=1,level=20)

   x,w = data.x, data.w
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(head[:wnames]))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

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

   title(head[:wnames][VarIndex_])

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
	plot(data, head, var; kwargs)

Wrapper over the plot function in matplotlib.
"""
function plot(data::Data, head::Dict, var::String; kwargs...)
   x,w = data.x, data.w
   VarIndex_ = findfirst(x->x==var,head[:wnames])

   c = plot(x, w[:,VarIndex_]; kwargs...)

   return c::Vector{PyCall.PyObject}
end

"""
	scatter(data, head, var; kwargs)

Wrapper over the scatter function in matplotlib.
"""
function scatter(data::Data, head::Dict, var::String; kwargs...)
   x,w = data.x, data.w
   VarIndex_ = findfirst(x->x==var,head[:wnames])

   c = scatter(x, w[:,VarIndex_]; kwargs...)

   return c::Vector{PyCall.PyObject}
end

"""
	contour(data, head, var, levels=0; plotrange, plotinterval, kwargs)

Wrapper over the contour function in matplotlib.
"""
function contour(data::Data, head::Dict, var::String, levels=0;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs=Dict())

   xi, yi, wi = getdata(data, head, var, plotrange, plotinterval)

   if levels != 0
      c = plt.contour(xi, yi, wi, levels; kwargs...)
   else
      c = plt.contour(xi, yi, wi; kwargs...)
   end

   return c::PyCall.PyObject
end

"""
	contourf(data, head, var, levels=0; plotrange, plotinterval, kwargs)

Wrapper over the contourf function in matplotlib.
"""
function contourf(data::Data, head::Dict, var::String, levels=0;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, kwargs=Dict())

   xi, yi, wi = getdata(data, head, var, plotrange, plotinterval)

   if levels != 0
      c = plt.contourf(xi, yi, wi, levels; kwargs...)
   else
      c = plt.contourf(xi, yi, wi; kwargs...)
   end

   return c::PyCall.PyObject
end

"""
	tricontourf(data, head, var; plotrange, plotinterval, kwargs)

Wrapper over the tricontourf function in matplotlib.
"""
function tricontourf(data::Data, head::Dict, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
   kwargs=Dict())

   x = data.x
   w = data.w

   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(head[:wnames]))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

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

   c = tricontourf(X, Y, W)
end

"""
	plot_trisurf(data::Data, head::Dict, var::String;
		plotrange::Vector{Float64}=[-Inf,Inf,-Inf,Inf], kwargs::Dict=Dict())

Wrapper over the plot_trisurf function in matplotlib.
"""
function plot_trisurf(data::Data, head::Dict, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], kwargs=Dict())

   x = data.x
   w = data.w

   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(head[:wnames]))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

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

   c = plot_trisurf(X, Y, W)
end


"""
	plot_surface(data, head, var; plotrange, plotinterval, kwargs)

Wrapper over the plot_surface function in matplotlib.
"""
function plot_surface(data::Data, head::Dict, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
   kwargs=Dict())

   xi, yi, wi = getdata(data, head, var, plotrange, plotinterval)

   c = plot_surface(xi, yi, wi; kwargs...)

   return c::PyCall.PyObject
end

"""
	streamplot(data, head, var; plotrange, plotinterval)

Wrapper over the streamplot function in matplotlib. Streamplot does not have
**kwargs in the API.
"""
function streamplot(data::Data, head::Dict, var::String;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1, density=1.0, color="")

   x, w = data.x, data.w
   VarStream  = split(var,";")
   wnames = lowercase.(head[:wnames])
   VarIndex1_ = findfirst(x->x==lowercase(VarStream[1]), wnames)
   VarIndex2_ = findfirst(x->x==lowercase(VarStream[2]), wnames)

   if head[:gencoord] # Generalized coordinates
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
	  Xi = [y for x in xi, y in yi]
	  Yi = [x for x in xi, y in yi]

	  interpolator = matplotlib.tri.LinearTriInterpolator(tr, w[:,1,VarIndex1_])
	  v1 = interpolator(Xi, Yi)

	  interpolator = matplotlib.tri.LinearTriInterpolator(tr, w[:,1,VarIndex2_])
	  v2 = interpolator(Xi, Yi)

   else # Cartesian coordinates
	  X = x[:,1,1]
	  Y = x[1,:,2]
	  if all(isinf.(plotrange))
		 Xi, Yi = X, Y
		 v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'
	  else
		 if plotrange[1] == -Inf plotrange[1] = minimum(X) end
		 if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
		 if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
		 if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

		 w1, w2 = w[:,:,VarIndex1_], w[:,:,VarIndex2_]

		 xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
		 yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

		 Xi = [i for i in xi, j in yi]
		 Yi = [j for i in xi, j in yi]

		 spline = Spline2D(X, Y, w1)
		 v1 = spline(Xi[:], Yi[:])
		 v1 = reshape(v1, size(Xi))'

		 spline = Spline2D(X, Y, w2)
		 v2 = spline(Xi[:], Yi[:])
		 v2 = reshape(v2, size(Xi))'
	  end
   end

   if isempty(color)
      c = streamplot(Xi, Yi, v1, v2; density=density)
   else
      c = streamplot(Xi, Yi, v1, v2; density=density, color=color)
   end
   return c::PyCall.PyObject
end

"""Prepare data for passing to plotting functions."""
function getdata(data, head, var, plotrange, plotinterval)
   x,w = data.x, data.w
   ndim = head[:ndim]

   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(head[:wnames]))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

   if head[:gencoord] # Generalized coordinates
      X = vec(x[:,:,1])
      Y = vec(x[:,:,2])
      W = vec(w[:,:,VarIndex_])

      if any(abs.(plotrange) .== Inf)
         if plotrange[1] == -Inf plotrange[1] = minimum(X) end
         if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
         if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
         if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
      end

      # Create grid values first.
      xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
      yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
      # Perform linear interpolation of the data (x,y) on grid(xi,yi)
      triang = matplotlib.tri.Triangulation(X,Y)
      interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
      Xi = [y for x in xi, y in yi]
	   Yi = [x for x in xi, y in yi]
      wi = interpolator(Xi, Yi)
   else # Cartesian coordinates
      if all(isinf.(plotrange))
         xi = x[:,:,1]
         yi = x[:,:,2]
         wi = w[:,:,VarIndex_]
      else
		   X = x[:,1,1]
		   Y = x[1,:,2]
         if plotrange[1] == -Inf plotrange[1] = minimum(X) end
         if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
         if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
         if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

         W = w[:,:,VarIndex_]

         xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
         yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

         spline = Spline2D(X, Y, W)
         Xi = [i for i in xi, j in yi]
         Yi = [j for i in xi, j in yi]
         wi = spline(Xi[:], Yi[:])
         wi = reshape(wi, size(Xi))'
      end
   end
   return xi, yi, wi
end

"""
	animatedata(filelist, func, (plotmode="contbar",
      plotrange=[-Inf Inf -Inf Inf],
      plotinterval=0.1))

Generate animations from data. This is basically calling plotdata function for
multiple snapshots. The main issue here is to determine the colorbar/axis range
in advance to avoid any jump in the movie.
"""
function animatedata(filelist::FileList,func::String;
   imin=1, imax=1, cut="", plotmode="contbar", plotrange=[-Inf,Inf,-Inf,Inf],
   plotinterval=0.1, verbose=true)
end


function animate(i,filelist)
   clf()
   fhead, d, flist = readdata(filelist.name,verbose=false,npict=i+1)
   plotdata(d[1],fhead[1],"p",plotmode="contbar")

   return gca()
end

"""
	convertVTK(head, data, connectivity, filename)

Convert 3D unstructured Tecplot data to VTK. Note that if using voxel type data
in VTK, the connectivity sequence is different from Tecplot.
"""
function convertVTK(head::Dict, data::Array{Float32,2},
   connectivity::Array{Int32,2}, filename="3DBATSRUS")

   nVar = length(head[:variables])

   points = @view data[1:head[:ndim],:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:ndim] == 3
      # PLT to VTK index_ = [1 2 4 3 5 6 8 7]
      for i = 1:2
         connectivity = swaprows(connectivity, 4*i-1, 4*i)
      end
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:ndim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid(filename, points, cells)

   for ivar = head[:ndim]+1:nVar
      if occursin("_x",head[:variables][ivar]) # vector
         var1 = @view data[ivar,:]
         var2 = @view data[ivar+1,:]
         var3 = @view data[ivar+2,:]
         namevar = replace(head[:variables][ivar], "_x"=>"")
         vtk_point_data(vtkfile, (var1, var2, var3), namevar)
      elseif occursin(r"(_y|_z)",head[:variables][ivar])
         continue
      else
         var = @view data[ivar,:]
         vtk_point_data(vtkfile, var, head[:variables][ivar])
      end
   end

   outfiles = vtk_save(vtkfile)
end

function get_var(data::Data, head::Dict, var::String)
   VarIndex_ = findfirst(x->x==var,head[:wnames])

   ndim = head[:ndim]
   if ndim == 1
      w = data.w[:,VarIndex_]
   elseif ndim == 2
      w = data.w[:,:,VarIndex_]
   elseif ndim == 3
      w = data.w[:,:,:,VarIndex_]
   end
   w
end

function get_vars(data::Data, head::Dict, Names::Vector{String})

   dict = Dict()
   for name in Names
      dict[name] = get_var(data, head, name)
   end

   Vars(dict)
end

Base.getproperty(p::Vars, name::Symbol) = getfield(p, :data)[String(name)]

function swaprows(X::Array{Int32,2}, i, j)
   m, n = size(X)
   if (1 <= i <= n) && (1 <= j <= n)
      for k = 1:n
        @inbounds X[i,k],X[j,k] = X[j,k],X[i,k]
      end
      return X
   else
      throw(BoundsError())
   end
end
