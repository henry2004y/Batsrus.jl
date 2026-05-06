module BatsrusMakieExt

using Batsrus, Printf
import Batsrus: findindex, hasunit, getunit, interp2d, plot_phase, plot_phase!,
    _resolve_alias
import Batsrus.UnitfulBatsrus
import Makie
using Makie: @L_str, heatmap!, poly!, Circle, Point2f, Figure, Axis, Colorbar,
    colgap!, DataAspect, save, delete!, Observable
using Unitful: ustrip
using DimensionalData: dims, val, X, Y, Z, lookup, rebuild
using PrecompileTools: @setup_workload, @compile_workload

include("typerecipe.jl")
include("makie_amrex.jl")

@setup_workload begin
    @compile_workload begin
        # Mock data for precompilation
        file = joinpath(@__DIR__, "../../test/precompile.out")
        bd = file |> load
        # Precompile common plotting functions
        Makie.heatmap(bd, "rho")
    end
end

end
