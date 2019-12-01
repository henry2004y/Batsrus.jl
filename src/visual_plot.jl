# Test of using user recipes from Plots.
#
# Hongyang Zhou, hyzhou@umich.edu 12/01/2019

using Plots

# Build a recipe which acts on a custom type.
# Notice that the function apply_recipe is returned.
# The recipe macro is just a convenience to build apply_recipe definitions.
@recipe function plot1d(data::Data, filehead::Dict, var::String)
   x,w = data.x, data.w
   VarIndex_ = findfirst(x->x==var,filehead[:wnames])
   y = w[:,VarIndex_]
   seriestype --> :path  # there is always an attribute dictionary `d` available...
   # If the user didn't specify a seriestype, we choose :path
   x, y
end

# that was pretty easy!
#plot(data[1], head[1], "p")
