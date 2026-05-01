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
    @compile_workload begin
        mktempdir() do tmpdir
            # Mock data for precompilation
            file = joinpath(@__DIR__, "../../test/precompile.out")
            animate([file]; var = :rho)
        end
    end
end

end
