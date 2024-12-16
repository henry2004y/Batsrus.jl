# Plot plasma β in a cut plane from 3D PIC outputs.
#
# Hongyang Zhou, hyzhou@umich.edu

using Batsrus, PyPlot, Printf

file = "3d_box.out"
data = load(file)

sequence = 65 # cut plane index from -
p_ = 18 # thermal pressure index

X, Y, Z = eachslice(data.x, dims=4)

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(3.3, 6)

W  = @view data.w[:,:,:,p_]
Bx = @view data.w[:,:,:,5]
By = @view data.w[:,:,:,6]
Bz = @view data.w[:,:,:,7]

cut1 = @view X[:,sequence,:]
cut2 = @view Z[:,sequence,:]
Pcut  = @view W[:,sequence,:]
Bxcut = @view Bx[:,sequence,:]
Bycut = @view By[:,sequence,:]
Bzcut = @view Bz[:,sequence,:]
PBcut = hypot.(Bx, By, Bz)

c = ax.contourf(cut1, cut2, Pcut./PBcut)
fig.colorbar(c; ax)
ax.axis("scaled")
title(data.head.wname[p_])

xlabel("x"); ylabel("z")

dim = [0.125, 0.013, 0.2, 0.045]
str = @sprintf "it=%d, time=%4.2f" data.head[:it] data.head[:time]
at = matplotlib.offsetbox.AnchoredText(str,
   loc="lower left", prop=Dict("size"=>8), frameon=true,
   bbox_to_anchor=(0., 1.),
   bbox_transform=ax.transAxes)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)
