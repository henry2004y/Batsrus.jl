# Makie recipe

using Batsrus
using Makie

include("utility.jl")

# 1D example
filename = "1d__raw_2_t25.60000_n00000258.out"
data = readdata(filename, dir="test/data")

var = "p"
VarIndex_ = findindex(data, var)

fig = Figure(resolution = (800, 600))
ax = Axis(fig[1,1], title="1D shocktube", xlabel="x", ylabel="Pressure [nPa]")
line1 = lines!(ax, data.x[:], data.w[:,VarIndex_], color = :red,
   label="Pressure")
axislegend()

fig
#=
# 2D example
#filename = "z=0_raw_1_t25.60000_n00000258.out"
#data = readdata(filename, dir="test/data")

var = "p"
plotrange = [-Inf,Inf,-Inf,Inf]
plotinterval = 0.1

fig = Figure(resolution = (800, 600))

xi, yi, wi = getdata(data, var, plotrange, plotinterval)

fig[1,1] = Axis(fig, title="2D shocktube")
c = heatmap!(fig[1,1], xi, yi, wi)
cbar = Colorbar(fig[1,2], c, label="Pressure [nPa]", width=20)

fig
=#