module BatsrusMakieUniformStreamlinesExt

using Batsrus, Printf
import Batsrus: hasunit, getunit, interp2d, animate
import Batsrus.UnitfulBatsrus
import Makie
using Makie: @L_str, heatmap!, poly!, Circle, Point2f, Figure, Axis, Colorbar,
    colgap!, DataAspect, save, delete!, Observable
using DimensionalData: dims, val
using UniformStreamlines: streamlines!, evenstream
using ProgressMeter

include("animate.jl")

end
