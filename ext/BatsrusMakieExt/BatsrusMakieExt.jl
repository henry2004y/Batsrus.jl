module BatsrusMakieExt

using Batsrus, Printf
import Batsrus: findindex, hasunit, getunit, interp2d, plot_phase, _resolve_alias
import Batsrus.UnitfulBatsrus
import Makie
using Makie: @L_str, heatmap!, poly!, streamlines!, Circle, Point2f, Figure, Axis, Colorbar,
    colgap!, DataAspect, save, delete!, Observable
using UniformStreamlines
using Unitful: ustrip

include("typerecipe.jl")
include("makie_amrex.jl")
include("animate.jl")

end
