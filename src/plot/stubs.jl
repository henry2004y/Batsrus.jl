export plotlogdata, plotgrid, cutplot, streamslice, plot_phase, plot_phase!, animate

"""
    plotlogdata(data, head::NamedTuple, func::AbstractString; plotmode = "line")

Plot log data. Stub to be implemented by plotting extensions.
"""
function plotlogdata end

"""
    plotgrid(bd::AbstractBatsrusData, ax=nothing; kwargs...)
    plotgrid(batl::Batl, ax=nothing; kwargs...)
    plotgrid(head, data, connectivity; ax=nothing, kwargs...)

Plot simulation mesh. Stub to be implemented by plotting extensions.
"""
function plotgrid end

"""
    cutplot(data, var, cut; kwargs...)

Plot a 2D slice of the data. Stub to be implemented by plotting extensions.
"""
function cutplot end

"""
    streamslice(data, var1, var2, cut; kwargs...)

Plot 2D streamlines. Stub to be implemented by plotting extensions.
"""
function streamslice end

"""
    plot_phase!(ax, data, x_variable, y_variable; kwargs...)

Plots the 2D phase space density for selected variables.
Stub to be implemented by plotting extensions.
"""
function plot_phase! end

"""
    animate(files::Vector{String}; kwargs...)

Save figures or create an animation from a sequence of SWMF output files.
Stub to be implemented by plotting extensions.

Notable behaviors:
* `vmin`/`vmax` set to `Inf` are determined by scanning **all** snapshots, so the
  color scale stays consistent across the animation.
* `plotinterval = nothing` auto-detects the finest grid resolution (see
  [`finest_resolution`](@ref)), which is essential for AMR/unstructured output.
* If `plot_kwargs` is not provided, the colormap is selected automatically:
  `RdBu_r` with symmetric limits for bipolar data, `turbo` otherwise.
"""
function animate end

"""
    plot_phase(data, args...; ax = nothing, kwargs...)

Plots the 2D phase space density for selected variables.
This is a wrapper around [`plot_phase!`](@ref) that allows passing the axis as a keyword
argument.
"""
function plot_phase(data, args...; ax = nothing, kwargs...)
    return plot_phase!(ax, data, args...; kwargs...)
end
