module BatsrusMakieExt

using Batsrus, Printf
import Batsrus: findindex, hasunit, getunit, interp2d
import Batsrus.UnitfulBatsrus
import Makie
using Makie: @L_str

import Batsrus: plot_phase, _resolve_alias

include("typerecipe.jl")
include("makie_amrex.jl")

end
