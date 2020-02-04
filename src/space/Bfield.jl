# Script for SWMF Galileo Flyby Simulation Comparison
#
# Hongyang Zhou, hyzhou@umich.edu 02/03/2020

using CSV, Dates, PyPlot, VisAna, Interpolations

#=
## Calculate differences & save to output file
Bdiff = Bobs - Bsim
# There`s a difference between 2 definitions of 2-norm.
norm2 = mean(abs(Bdiff(~isnan(Bdiff))).^2)^(1/2)
# Find the index and value of the best parameter set
[val, idx] = min(norm2)
=#


# Parameters
flyby =  8   # [1,2,7,8,28,29]
DoPlot = true  # Plot output
DoSave = false  # Save norm2 number
firstpict = 1 # first snapshot to pick
lastpict  = 1 # last snapshot to pick
fileGathered = false # Input data in 1 file or multiple files

# Read observation data
dir = "/Users/hyzhou/Documents/research/Ganymede/Galileo"
flybyfile = "Galileo_G"*string(flyby)*"_flyby_MAG.dat"

df = CSV.File(joinpath(dir,flybyfile);header=2,delim=" ",ignorerepeated=true)
t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
   floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )
BobsStrength = @. sqrt(df.Bx^2 + df.By^2 + df.Bz^2)

# Read simulation data
if flyby == 8
   filename=""
   # Select the starting time for synthetic satellite
   tStart = DateTime(1997,5,7,15,48,7)
elseif flyby == 28
   filename=""
   # Select the starting time for synthetic satellite
   tStart = DateTime(2000,5,20,10,0,12)
end

head, data, list = readdata(filename,
   dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia",
   npict=npict)

# Interpolate simulation data to observation data
nx, nw = head[1][:nx], head[1][:nw]
npict = list[1][:npictinfiles]

if firstpict > npict
   @error "firstpict out of range!"
end

npict += 1 - firstpict

kStart = findfirst(tStart .- t .< Second(0))

tEnd  = tStart + Second(npict-1)
kEnd = findfirst(tEnd .- t .< Second(0))

nSim = kStart + (length(t) - kEnd + 1) + npict - 2
BxSim = Vector{Float32}(undef,nSim)
BySim = Vector{Float32}(undef,nSim)
BzSim = Vector{Float32}(undef,nSim)
tSim = Vector{DateTime}(undef,nSim)
tSim[1:kStart] = t[1:kStart]

function mylinspace(d1::DateTime, d2::DateTime; n=2)
    Δ = d2 - d1
    T = typeof(Δ)
    δ = T(round(Int, Dates.value(Δ)/(n - 1)))
    d2 = d1 + δ*(n - 1)
    return d1:δ:d2
end

# This is not that accurate, but just an approximation
tSim[kStart:kStart+npict-1] = mylinspace(tStart, tEnd, n=npict)
tSim[kStart+npict-1:end] = t[kEnd:end]

# Extract the magnetic field along the trajectory before the starting time
# in one snapshot
head, data = read_data(filename, npict=firstpict)

# Interpolate simulation data to observation data
x = data[1].x[:,:,:,1]
y = data[1].x[:,:,:,2]
z = data[1].x[:,:,:,3]

bx_ = findfirst(x->x=="bx", head[1][:wnames])
by_ = findfirst(x->x=="by", head[1][:wnames])
bz_ = findfirst(x->x=="bz", head[1][:wnames])
Bx = @view data[1].w[:,:,:,bx_]
By = @view data[1].w[:,:,:,by_]
Bz = @view data[1].w[:,:,:,bz_]

# Gridded interpolation
knots = (x, y, z)
itpBx = interpolate(knots, Bx, Gridded(Linear()))
BxSim[1:kStart] = itpBx.(df.X[1:kStart],df.Y[1:kStart],df.Z[1:kStart])

itpBy = interpolate(knots, By, Gridded(Linear()))
BySim[1:kStart] = itpBy.(df.X[1:kStart],df.Y[1:kStart],df.Z[1:kStart])

itpBz = interpolate(knots, Bz, Gridded(Linear()))
BzSim[1:kStart] = itpBz.(df.X[1:kStart],df.Y[1:kStart],df.Z[1:kStart])



# Extract the magnetic field along the trajectory from multiple snapshot
for ipict = 1:npict-1
   println("ipict=$(round(ipict,2))",)
   head, data = readdata(filename, npict=firstpict+ipict)
   # Interpolate simulation data to observation data
   x = data[1].x[:,:,:,1]
   y = data[1].x[:,:,:,2]
   z = data[1].x[:,:,:,3]

   Bx = @view data[1].w[:,:,:,bx_]
   By = @view data[1].w[:,:,:,by_]
   Bz = @view data[1].w[:,:,:,bz_]

   k = argmin( abs(timesim(kStart+ipict) - t) )

   itpBx = interpolate(knots, Bx, Gridded(Linear()))
   BxSim[kStart+ipict] = itpBx(df.X[k],df.Y[k],df.Z[k])

   itpBy = interpolate(knots, By, Gridded(Linear()))
   BySim[kStart+ipict] = itpBy(df.X[k],df.Y[k],df.Z[k])

   itpBz = interpolate(knots, Bz, Gridded(Linear()))
   BzSim[kStart+ipict] = itpBz(df.X[k],df.Y[k],df.Z[k])
end


BxSim[kStart+npict:end,1] = itpBx.(df.X[kEnd+1:end],df.Y[kEnd+1:end],df.Z[kEnd+1:end])
BySim[kStart+npict:end,2] = itpBy.(df.X[kEnd+1:end],df.Y[kEnd+1:end],df.Z[kEnd+1:end])
BzSim[kStart+npict:end,3] = itpBz.(df.X[kEnd+1:end],df.Y[kEnd+1:end],df.Z[kEnd+1:end])

BsimStrength = @. sqrt(BxSim^2 + BySim^2 + BzSim^2)

# Visualization
figure(figsize=(12,5))
formatter = matplotlib.dates.DateFormatter("%H:%M")
ax1 = subplot(411)
plot(t, df.Bx, label=L"B_x")
plot(t, BxSim)
ax1.xaxis.set_major_formatter(formatter)
legend()
ax2 = subplot(412)
plot(t, df.By, label=L"B_y")
plot(t, BySim)
ax2.xaxis.set_major_formatter(formatter)
legend()
ax3 = subplot(413)
plot(t, df.Bz, label=L"B_z")
plot(t, BzSim)
ax3.xaxis.set_major_formatter(formatter)
legend()
ax4 = subplot(414)
plot(t, BobsStrength, label=L"B")
plot(t, BsimStrength)
ax4.xaxis.set_major_formatter(formatter)
legend()
tight_layout()

# Save simulation line data
if DoSave

end



"""
	plotBSteady(flyby=8, DoSave=true)

Plot quasi-steady state magnetic field comparison. Second norm is calculated to
quantitatively estimate the difference.
`flyby` is an integer in [1,2,7,8,28,29].
`DoSave` decides
"""
function plotBSteady(flyby=8, DoSave=false)
   # Read observation data
   dir = "/Users/hyzhou/Documents/research/Ganymede/Galileo"
   flybyfile = "Galileo_G"*string(flyby)*"_flyby_MAG.dat"

   df = CSV.File(joinpath(dir,flybyfile);header=2,delim=" ",ignorerepeated=true)
   t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
   	floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )
   BobsStrength = @. sqrt(df.Bx^2 + df.By^2 + df.Bz^2)

   # Read simulation data
   filename="box_var_4_n00080000.out"
   npict = 1 # Remember to change this for different runs!
   head, data = readdata(filename,
      dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia",
      npict=npict)

   # Interpolate simulation data to observation data
   nx, nw = head[1][:nx], head[1][:nw]

   x = data[1].x[:,1,1,1]
   y = data[1].x[1,:,1,2]
   z = data[1].x[1,1,:,3]

   w = data[1].w
   # Find the correct index of variables
   bx_ = findfirst(x->x=="bx", head[1][:wnames])
   by_ = findfirst(x->x=="by", head[1][:wnames])
   bz_ = findfirst(x->x=="bz", head[1][:wnames])
   Bx = @view data[1].w[:,:,:,bx_]
   By = @view data[1].w[:,:,:,by_]
   Bz = @view data[1].w[:,:,:,bz_]

   # Gridded interpolation
   knots = (x, y, z)
   itpBx = interpolate(knots, Bx, Gridded(Linear()))
   BxSim = itpBx.(df.X,df.Y,df.Z)

   itpBy = interpolate(knots, By, Gridded(Linear()))
   BySim = itpBy.(df.X,df.Y,df.Z)

   itpBz = interpolate(knots, Bz, Gridded(Linear()))
   BzSim = itpBz.(df.X,df.Y,df.Z)

   BsimStrength = @. sqrt(BxSim^2 + BySim^2 + BzSim^2)

   # Visualization
   figure(figsize=(12,5))
   formatter = matplotlib.dates.DateFormatter("%H:%M")
   ax1 = subplot(411)
   plot(t, df.Bx, label=L"B_x")
   plot(t, BxSim)
   ax1.xaxis.set_major_formatter(formatter)
   legend()
   ax2 = subplot(412)
   plot(t, df.By, label=L"B_y")
   plot(t, BySim)
   ax2.xaxis.set_major_formatter(formatter)
   legend()
   ax3 = subplot(413)
   plot(t, df.Bz, label=L"B_z")
   plot(t, BzSim)
   ax3.xaxis.set_major_formatter(formatter)
   legend()
   ax4 = subplot(414)
   plot(t, BobsStrength, label=L"B")
   plot(t, BsimStrength)
   ax4.xaxis.set_major_formatter(formatter)
   legend()
   tight_layout()

   if DoSave # Save the plot

   end

end

"""


Plot time-accurate magnetic field comparisons.
This function is able to handle any npict number of snapshots. If the simulation
time interval is smaller than observation (which is usually the case), than the
magnetic field along the trajectoryy before time_start is obtained from first
single snapshot ipict=1, and B after time_end is obtained from the last single
snapshot ipict=npict. In this way, I avoid the jump due to periodic
interpolation.
`flyby` and `time_start` need to be modified when switching between flybys.
"""
function plotBTimeAccurate()

end
