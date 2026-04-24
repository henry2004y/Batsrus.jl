export plotlogdata, cutplot, streamslice, plot_phase, plot_phase!, animate

function plotlogdata end
function cutplot end
function streamslice end
function plot_phase! end
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
