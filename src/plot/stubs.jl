export plotlogdata, cutplot, streamslice, plot_phase

function plotlogdata end
function cutplot end
function streamslice end
function plot_phase! end

function plot_phase(data, args...; ax = nothing, kwargs...)
   plot_phase!(ax, data, args...; kwargs...)
end
