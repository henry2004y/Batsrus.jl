# Plotting functionalities.

using PyPlot

export plotlogdata, plot, scatter, contour, contourf, plot_surface, tripcolor,
    tricontourf, plot_trisurf, streamplot, streamslice, quiver, cutplot, pcolormesh

"""
    plotlogdata(data, head, func; plotmode="line")

Plot information from log file.

# Input arguments

  - `data::Array`: output data.
  - `head::NamedTuple`: header info.
  - `func::String`: variables for plotting.
  - `plotmode::String`: type of plotting ["line","scatter"].
"""
function plotlogdata(data, head::NamedTuple, func::AbstractString; plotmode = "line")
    vars = split(func)
    plotmode = split(plotmode)

    for (ivar, var) in enumerate(vars)
        varIndex_ = findfirst(x -> lowercase(x) == lowercase(var), head.variable)
        isnothing(varIndex_) && error("$(var) not found in file header variables!")

        figure()
        if plotmode[ivar] == "line"
            plot(data[1, :], data[varIndex_, :])
        elseif plotmode[ivar] == "scatter"
            scatter(data[1, :], data[varIndex_, :])
        else
            throw(ArgumentError("unknown plot mode $(plotmode[ivar])!"))
        end
        xlabel(head.variable[1])
        ylabel(head.variable[varIndex_])
        title("log file data")
    end
    return
end

"""
    plotgrid(bd::BATS, var, ax=nothing; kwargs...)

Plot 2D mesh.
"""
function plotgrid(
        bd::BatsrusIDL{2, TV},
        func::AbstractString,
        ax = nothing;
        kwargs...
    ) where {TV}
    if isnothing(ax)
        ax = plt.gca()
    end

    # This does not take subdomain plot into account!
    X, Y = eachslice(bd.x, dims = 3)
    scatter(X, Y, marker = ".", alpha = 0.6)
    title("Grid illustration")

    xlabel(bd.head.wname[1])
    ylabel(bd.head.wname[2])
    add_time_iteration!(bd, ax)

    return
end

"""
    cutplot(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

2D plane cut pcolormesh of 3D box data. `sequence` is the index along `dir`.
"""
function cutplot(
        bd::BatsrusIDL{3, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], dir = "x", sequence = 1
    ) where {TV}
    x, w = bd.x, bd.w
    varIndex_ = findindex(bd, var)

    X, Y, Z = eachslice(x, dims = 4)

    W = @view w[:, :, :, varIndex_]

    if dir == "x"
        cut1 = @view X[sequence, :, :]
        cut2 = @view Y[sequence, :, :]
        W = @view W[sequence, :, :]
    elseif dir == "y"
        cut1 = @view X[:, sequence, :]
        cut2 = @view Z[:, sequence, :]
        W = @view W[:, sequence, :]
    elseif dir == "z"
        cut1 = @view X[:, :, sequence]
        cut2 = @view Y[:, :, sequence]
        W = @view W[:, :, sequence]
    end

    if !all(isinf.(plotrange))
        cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
    end
    if isnothing(ax)
        ax = plt.gca()
    end
    c = ax.pcolormesh(cut1, cut2, W)

    title(bd.head.wname[varIndex_])

    if dir == "x"
        xlabel("y")
        ylabel("z")
    elseif dir == "y"
        xlabel("x")
        ylabel("z")
    elseif dir == "z"
        xlabel("x")
        ylabel("y")
    end

    add_time_iteration!(bd, ax)

    return c
end

"""
    streamslice(data::BATS, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], dir="x",
    	 sequence=1; kwargs...)

Plot streamlines on 2D slices of 3D box data. Variable names in `var` string must be
separated with `;`.
"""
function streamslice(
        bd::BatsrusIDL{3, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], dir = "x", sequence = 1, kwargs...
    ) where {TV}
    x, w = bd.x, bd.w
    varstream = split(var, ";")
    var1_ = findindex(bd, varstream[1])
    var2_ = findindex(bd, varstream[2])

    X, Y, Z = eachslice(x, dims = 4)

    v1 = @view w[:, :, :, var1_]
    v2 = @view w[:, :, :, var2_]

    if dir == "x"
        cut1 = @view X[sequence, :, :]
        cut2 = @view Y[sequence, :, :]
        v1 = v1[sequence, :, :]
        v2 = v2[sequence, :, :]
    elseif dir == "y"
        cut1 = @view X[:, sequence, :]
        cut2 = @view Z[:, sequence, :]
        v1 = v1[:, sequence, :]
        v2 = v2[:, sequence, :]
    elseif dir == "z"
        cut1 = @view X[:, :, sequence]
        cut2 = @view Y[:, :, sequence]
        v1 = v1[:, :, sequence]
        v2 = v2[:, :, sequence]
    end

    if !all(isinf.(plotrange))
        cut1, cut2, v1, v2 = subsurface(cut1, cut2, v1, v2, plotrange)
    end

    xi = range(
        cut1[1, 1], stop = cut1[end, 1],
        step = (cut1[end, 1] - cut1[1, 1]) / (size(cut1, 1) - 1)
    )
    yi = range(
        cut2[1, 1], stop = cut2[1, end],
        step = (cut2[1, end] - cut2[1, 1]) / (size(cut2, 2) - 1)
    )

    if isnothing(ax)
        ax = plt.gca()
    end
    s = ax.streamplot(xi, yi, v1', v2'; kwargs...)

    if dir == "x"
        xlabel("y")
        ylabel("z")
    elseif dir == "y"
        xlabel("x")
        ylabel("z")
    elseif dir == "z"
        xlabel("x")
        ylabel("y")
    end

    return s
end

"""
    plot(bd::BATS{1, TV, T}, var, ax=nothing; kwargs...)

Wrapper over `plot` in matplotlib. Plot 1D outputs.
"""
function PyPlot.plot(
        bd::BatsrusIDL{1, TV},
        var::AbstractString,
        ax = nothing;
        kwargs...
    ) where {TV}
    x, w = bd.x, bd.w
    varIndex_ = findindex(bd, var)
    if isnothing(ax)
        ax = plt.gca()
    end

    c = ax.plot(x, w[:, varIndex_]; kwargs...)

    xlabel("x")
    ylabel("$(var)")
    add_time_iteration!(bd, ax)

    return c
end

"""
    scatter(data, var, ax=nothing; kwargs...)

Wrapper over `scatter` in matplotlib.
"""
function PyPlot.scatter(
        bd::BatsrusIDL{1, TV},
        var::AbstractString,
        ax = nothing;
        kwargs...
    ) where {TV}
    x, w = bd.x, bd.w
    varIndex_ = findindex(bd, var)
    if isnothing(ax)
        ax = plt.gca()
    end

    return c = ax.scatter(x, w[:, varIndex_]; kwargs...)
end

"""
    contour(data, var, ax=nothing; levels=0, plotrange=[-Inf,Inf,-Inf,Inf],
    	 plotinterval=0.1, innermask=false, kwargs...)

Wrapper over `contour` in matplotlib.
"""
function PyPlot.contour(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        levels = 0,
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1, innermask = false, rbody = 1.0,
        kwargs...
    ) where {TV}
    Xi, Yi, Wi = interp2d(bd, var, plotrange, plotinterval; innermask, rbody)
    if isnothing(ax)
        ax = plt.gca()
    end

    if levels != 0
        c = ax.contour(Xi, Yi, Wi, levels; kwargs...)
    else
        c = ax.contour(Xi, Yi, Wi; kwargs...)
    end
    add_titles!(bd, var, ax)

    return c
end

"""
    contourf(data, var, ax=nothing; levels=0, plotrange=[-Inf,Inf,-Inf,Inf],
    	 plotinterval=0.1, innermask=false, add_colorbar=true, vmin=-Inf, vmax=Inf,
    	 colorscale=:linear, kwargs...)

Wrapper over `contourf` in matplotlib. See [`interp2d`](@ref) for some related keywords.
"""
function PyPlot.contourf(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        levels::Int = 0,
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1, innermask = false, rbody = 1.0,
        add_colorbar = true, vmin = -Inf, vmax = Inf, colorscale = :linear, kwargs...
    ) where {TV}
    Xi, Yi, Wi = interp2d(bd, var, plotrange, plotinterval; innermask, rbody)
    if isnothing(ax)
        ax = plt.gca()
    end

    norm = set_colorbar(colorscale, vmin, vmax, Wi)
    if levels != 0
        c = ax.contourf(Xi, Yi, Wi, levels; norm, kwargs...)
    else
        c = ax.contourf(Xi, Yi, Wi; norm, kwargs...)
    end
    add_colorbar && colorbar(c; ax, fraction = 0.04, pad = 0.02)
    add_titles!(bd, var, ax)

    return c
end

"""
     tricontourf(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], kwargs...)

Wrapper over `tricontourf` in matplotlib.
"""
function PyPlot.tricontourf(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], kwargs...
    ) where {TV}
    x, w = bd.x, bd.w
    varIndex_ = findindex(bd, var)

    X = vec(x[:, :, 1])
    Y = vec(x[:, :, 2])
    W = vec(w[:, :, varIndex_])

    #TODO This needs improvement.
    if !all(isinf.(plotrange))
        xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
            Y .> plotrange[3] .& Y .< plotrange[4]
        X = X[xyIndex]
        Y = Y[xyIndex]
        W = W[xyIndex]
    end
    if isnothing(ax)
        ax = plt.gca()
    end

    c = ax.tricontourf(X, Y, W; kwargs...)

    add_titles!(bd, var, ax)

    return c
end

function PyPlot.triplot(
        bd::BatsrusIDL{2, TV}, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf],
        kwargs...
    ) where {TV}
    X = vec(bd.x[:, :, 1])
    Y = vec(bd.x[:, :, 2])
    triang = PyPlot.matplotlib.tri.Triangulation(X, Y)
    #TODO This needs improvement.
    if !all(isinf.(plotrange))
        xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
            Y .> plotrange[3] .& Y .< plotrange[4]
        X = X[xyIndex]
        Y = Y[xyIndex]
    end
    if isnothing(ax)
        ax = plt.gca()
    end

    return c = ax.triplot(triang, kwargs...)
end

"""
    plot_trisurf(data::BATS, var::String, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf],
    	 kwargs...)

Wrapper over `plot_trisurf` in matplotlib.
"""
function PyPlot.plot_trisurf(
        bd::BatsrusIDL{2, TV}, var::AbstractString;
        plotrange = [-Inf, Inf, -Inf, Inf], kwargs...
    ) where {TV}
    x, w = bd.x, bd.w
    varIndex_ = findindex(bd, var)

    X = vec(x[:, :, 1])
    Y = vec(x[:, :, 2])
    W = vec(w[:, :, varIndex_])

    #TODO This needs improvement.
    if !all(isinf.(plotrange))
        xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
            Y .> plotrange[3] .& Y .< plotrange[4]
        X = X[xyIndex]
        Y = Y[xyIndex]
        W = W[xyIndex]
    end
    if isnothing(ax)
        ax = plt.gca()
    end

    ax = plt.figure().add_subplot(projection = "3d")
    c = ax.plot_trisurf(X, Y, W)

    add_titles!(bd, var, ax)

    return c
end

"""
    plot_surface(data, var; plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
    	 innermask=false, kwargs...)

Wrapper over `plot_surface` in matplotlib.
"""
function PyPlot.plot_surface(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1, innermask = false, rbody = 1.0,
        kwargs...
    ) where {TV}
    if isnothing(ax)
        ax = plt.gca()
    end
    xi, yi, Wi = interp2d(bd, var, plotrange, plotinterval; innermask, rbody)
    Xi, Yi = meshgrid(xi, yi)

    c = plot_surface(Xi, Yi, Wi; kwargs...)

    add_titles!(bd, var, ax)

    return c
end

"""
    pcolormesh(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf],
    	 plotinterval=0.1, innermask=false, vmin=-Inf, vmax=Inf, colorscale=:linear,
    	 add_colorbar=true, kwargs...)

# Keywords

  - `rbody=1.0`: inner body radius.
  - `colorscale::Symbol`: colormap scale from [`:linear`, `:log`].
  - `add_colorbar=true`: turn on colorbar.

Wrapper over `pcolormesh` in matplotlib.
"""
function PyPlot.pcolormesh(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1, innermask = false, rbody = 1.0,
        vmin = -Inf, vmax = Inf, colorscale = :linear, add_colorbar = true, kwargs...
    ) where {TV}
    xi, yi, Wi = interp2d(bd, var, plotrange, plotinterval; innermask, rbody)

    if isnothing(ax)
        ax = plt.gca()
    end

    norm = set_colorbar(colorscale, vmin, vmax, Wi)
    c = ax.pcolormesh(xi, yi, Wi; norm, kwargs...)

    add_colorbar && colorbar(c; ax, fraction = 0.04, pad = 0.02)

    add_titles!(bd, var, ax)

    return c
end

"""
    tripcolor(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf],
    	 plotinterval=0.1, innermask=false, kwargs...)

Wrapper over `tripcolor` in matplotlib.
"""
function PyPlot.tripcolor(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], innermask = false, kwargs...
    ) where {TV}
    x, w = bd.x, bd.w

    varIndex_ = findindex(bd, var)

    X, Y = eachslice(x, dims = 3)
    adjust_plotrange!(plotrange, extrema(X), extrema(Y))
    W = vec(w[:, :, varIndex_])

    adjust_plotrange!(plotrange, extrema(X), extrema(Y))
    triang = PyPlot.matplotlib.tri.Triangulation(vec(X), vec(Y))

    # Mask off unwanted triangles at the inner boundary.
    if innermask
        varIndex_ = findlast(x -> x == "rbody", bd.head.param)
        isnothing(varIndex_) && error("rbody not found in file header parameters!")
        ParamIndex_ = varIndex_ - 2 - bd.head.nw
        r2 = bd.head.eqpar[ParamIndex_]^2

        ids = triang.triangles .+ Int32(1)
        mask = Vector{Bool}(undef, size(ids, 1))
        for i in axes(ids, 1)
            xmean = sum(@view X[ids[i, :]]) / 3
            ymean = sum(@view Y[ids[i, :]]) / 3
            xmean^2 + ymean^2 < r2 && (mask[i] = true)
        end

        triang.set_mask(mask)
    end
    if isnothing(ax)
        ax = plt.gca()
    end

    c = ax.tripcolor(triang, W; kwargs...)

    ax.set_xlim(plotrange[1], plotrange[2])
    ax.set_ylim(plotrange[3], plotrange[4])

    add_titles!(bd, var, ax)

    return c
end

"""
    imshow(data, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
    	 innermask=false, colorscale=:linear, add_colorbar=true, kwargs...)

# Keywords

  - `rbody=1.0`: inner body radius.
  - `colorscale::Symbol`: colormap scale from [`:linear`, `:log`].
  - `vmin`: minimum value of colorbar.
  - `vmax`: maximum value of colorbar.
  - `add_colorbar=true`: turn on colorbar.

Wrapper over `imshow` in matplotlib. For large matrices, this is faster than `pcolormesh`.
"""
function PyPlot.imshow(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1, innermask = false, rbody = 1.0,
        add_colorbar = true, vmin = -Inf, vmax = Inf, colorscale = :linear, kwargs...
    ) where {
        TV,
    }
    xi, yi, Wi = interp2d(bd, var, plotrange, plotinterval; innermask, rbody)

    if isnothing(ax)
        ax = plt.gca()
    end

    norm = set_colorbar(colorscale, vmin, vmax, Wi)
    c = ax.imshow(
        Wi; extent = [xi[1], xi[end], yi[1], yi[end]],
        origin = "lower", aspect = "auto", interpolation = "nearest", norm, kwargs...
    )

    add_colorbar && colorbar(c; ax, fraction = 0.04, pad = 0.02)

    add_titles!(bd, var, ax)

    return c
end

"""
    streamplot(bd, var, ax=nothing; plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1,
    	 kwargs...)

Wrapper over `streamplot` in matplotlib.
"""
function PyPlot.streamplot(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = Inf, kwargs...
    ) where {TV}
    xi, yi, v1, v2 = _getvector(bd, var; plotrange, plotinterval)

    if isnothing(ax)
        ax = plt.gca()
    end

    return ax.streamplot(xi, yi, v1, v2; kwargs...)
end

function _getvector(
        bd::BatsrusIDL{2, TV}, var::AbstractString;
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = Inf
    ) where {TV}
    x, w = bd.x, bd.w
    varstream = split(var, ";")
    var1_ = findfirst(x -> lowercase(x) == lowercase(varstream[1]), bd.head.wname)
    var2_ = findfirst(x -> lowercase(x) == lowercase(varstream[2]), bd.head.wname)
    if isinf(plotinterval)
        plotinterval = (x[end, 1, 1] - x[1, 1, 1]) / size(x, 1)
    end
    if bd.head.gencoord # generalized coordinates
        X, Y = vec(x[:, :, 1]), vec(x[:, :, 2])
        adjust_plotrange!(plotrange, extrema(X), extrema(Y))

        # Create grid values first.
        xi = range(Float64(plotrange[1]), stop = Float64(plotrange[2]), step = plotinterval)
        yi = range(Float64(plotrange[3]), stop = Float64(plotrange[4]), step = plotinterval)

        # Is there a triangulation method in Julia?
        tr = PyPlot.matplotlib.tri.Triangulation(X, Y)
        Xi, Yi = meshgrid(xi, yi)

        interpolator = PyPlot.matplotlib.tri.LinearTriInterpolator(tr, w[:, 1, var1_])
        v1 = interpolator(Xi, Yi)

        interpolator = PyPlot.matplotlib.tri.LinearTriInterpolator(tr, w[:, 1, var2_])
        v2 = interpolator(Xi, Yi)
    else # Cartesian coordinates
        # Convert to Float64 to satisfy the equal space checking in streamplot.py
        xrange, yrange = get_range(bd)
        if all(isinf.(plotrange))
            xi = xrange
            yi = yrange
            v1 = w[:, :, var1_].data'
            v2 = w[:, :, var2_].data'
        else
            adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

            w1, w2 = @views w[:, :, var1_], w[:, :, var2_]
            xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
            yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
            xyrange = (xrange, yrange)
            interp1 = cubic_spline_interpolation(xyrange, w1)
            v1 = [interp1(i, j) for j in yi, i in xi]

            interp2 = cubic_spline_interpolation(xyrange, w2)
            v2 = [interp2(i, j) for j in yi, i in xi]
        end
    end

    return xi, yi, v1, v2
end

"""
    quiver(data, var, ax=nothing; stride=10, kwargs...)

Wrapper over `quiver` in matplotlib. Only supports Cartesian grid for now.
"""
function PyPlot.quiver(
        bd::BatsrusIDL{2, TV}, var::AbstractString, ax = nothing;
        stride::Integer = 10, kwargs...
    ) where {TV}
    x, w = bd.x, bd.w
    VarQuiver = split(var, ";")
    var1_ = findindex(bd, VarQuiver[1])
    var2_ = findindex(bd, VarQuiver[2])

    @views X, Y = x[:, 1, 1], x[1, :, 2]
    @views v1, v2 = w[:, :, var1_]', w[:, :, var2_]'

    @views Xq, Yq = X[1:stride:end], Y[1:stride:end]
    v1q = @view v1[1:stride:end, 1:stride:end]
    v2q = @view v2[1:stride:end, 1:stride:end]

    if isnothing(ax)
        ax = plt.gca()
    end

    return ax.quiver(Xq, Yq, v1q, v2q; kwargs...)
end

"""
Set colorbar norm and ticks.
"""
function set_colorbar(colorscale, vmin, vmax, data = [1.0])
    if colorscale == :linear || any(<(0), data)
        colorscale == :log && @warn "Nonpositive data detected: use linear scale instead!"
        vmin = isinf(vmin) ? minimum(x -> isnan(x) ? +Inf : x, data) : vmin
        vmax = isinf(vmax) ? maximum(x -> isnan(x) ? -Inf : x, data) : vmax
        cnorm = PyPlot.matplotlib.colors.Normalize(vmin, vmax)
    elseif colorscale == :symlog
        cnorm = PyPlot.matplotlib.colors.SymLogNorm(
            linthresh = 0.03, linscale = 0.75; vmin, vmax
        )
    else # logarithmic
        datapositive = data[data .> 0.0]
        vmin = isinf(vmin) ? minimum(datapositive) : vmin
        vmax = isinf(vmax) ? maximum(x -> isnan(x) ? -Inf : x, data) : vmax
        cnorm = PyPlot.matplotlib.colors.LogNorm(vmin, vmax)
    end

    return cnorm
end

function add_titles!(bd::BatsrusIDL, var, ax)
    varIndex_ = findindex(bd, var)
    title(bd.head.wname[varIndex_])

    xlabel(bd.head.coord[1])
    ylabel(bd.head.coord[2])
    return add_time_iteration!(bd, ax)
end

function add_time_iteration!(bd::BatsrusIDL, ax)
    str = @sprintf "it=%d, time=%4.2f" bd.head.it bd.head.time
    at = PyPlot.matplotlib.offsetbox.AnchoredText(
        str,
        loc = "lower left", prop = Dict("size" => 8), frameon = true,
        bbox_to_anchor = (0.0, 1.0), bbox_transform = ax.transAxes
    )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    return ax.add_artist(at)
end

include("pyplot_amrex.jl")
