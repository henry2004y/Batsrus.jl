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
using PrecompileTools: @setup_workload, @compile_workload

include("animate.jl")

@setup_workload begin
    # Mock data for precompilation
    file = joinpath(@__DIR__, "../../test/precompile.out")
    @compile_workload begin
        # Precompile common plotting functions
        animate([file], "rho"; showplot = false, save_as = "")
    end
end

end
