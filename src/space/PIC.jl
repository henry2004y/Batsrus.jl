# Plotting PIC variables from box outputs.
#
# E: [μV/m]
# B: [nT]
# V: [km/s]
# ρ: [amu/cc]
# P: [nPa]
#
# Hongyang Zhou, hyzhou@umich.edu 02/06/2020

using VisAna, PyPlot

# For precise colorbar control
using PyCall
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
inset_axes = inset_locator.inset_axes
# For specifying the zero point in the colorbar
DN = matplotlib.colors.DivergingNorm

#dir = "/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia"
dir = "/Users/hyzhou"
#fnameField = "3d_var_region0_0_t00001640_n00020369.out"
#fnameField = "3d_var_region0_0_t00001520_n00004093.out"
#fnameField = "3d_var_region0_0_t00001523_n00004732.out"
#fnameField = "3d_var_region0_0_t00001524_n00004933.out"
#fnameField = "3d_var_region0_0_t00001525_n00005126.out"
fnameField = "3d_var_region0_0_t00001527_n00005528.out"
#fnameField = "3d_var_region0_0_t00001740_n00032921.out"
#fnameField = "3d_var_region0_0_t00001750_n00035011.out"
#fnameField = "3d_var_region0_0_t00001800_n00037224.out"
#fnameField = "3d_var_region0_0_t00001810_n00039298.out"
#fnameField = "3d_var_region0_0_t00001820_n00041392.out"
#fnameField = "3d_var_region0_0_t00001830_n00043324.out"
#fnameField = "3d_var_region0_0_t00001840_n00045319.out"
#fnameField = "3d_var_region0_0_t00001850_n00047357.out"
#fnameField = "3d_var_region0_0_t00001900_n00049293.out"
#fnameField = "3d_var_region0_0_t00001910_n00051301.out"

head, data = readdata(fnameField, dir=dir)

me = head[1][:eqpar][1]
qe = head[1][:eqpar][2]
mi = head[1][:eqpar][3]
qi = head[1][:eqpar][4]
const kB = 1.38064852e-23 # [m^2 kg s^-2 K^-1]
const vAlfven = 253. # reference Alfven velocity, [km/s]
const B₀ = √((-10.)^2+(-6.)^2+(-86.)^2)
const E₀ = vAlfven*B₀ # [μV/m]
const ρ₀ = 56.0     # [amu/cc]
const J₀ = 4.0*vAlfven
#const T₀ = vAlfven^2
const T₀ = 0.2/4 # Pe/n₀

#plotrange = [-2.05, -1.75, -0.5, 0.5]
plotrange = [-2.12, -1.75, -0.65, 0.6]
#plotrange=[-Inf, Inf, -Inf, Inf]
cI = 129 # plane cut index


X, Z, ρe = cutdata(data[1],head[1],"rhoS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, ρi = cutdata(data[1],head[1],"rhoS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Bx = cutdata(data[1],head[1],"Bx",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, By = cutdata(data[1],head[1],"By",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Bz = cutdata(data[1],head[1],"Bz",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Ex = cutdata(data[1],head[1],"Ex",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Ey = cutdata(data[1],head[1],"Ey",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Ez = cutdata(data[1],head[1],"Ez",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Uxe= cutdata(data[1],head[1],"uxS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Uye= cutdata(data[1],head[1],"uyS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Uze= cutdata(data[1],head[1],"uzS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Uxi= cutdata(data[1],head[1],"uxS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Uyi= cutdata(data[1],head[1],"uyS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Uzi= cutdata(data[1],head[1],"uzS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)

X, Z, Pxxe = cutdata(data[1],head[1],"pXXS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pyye = cutdata(data[1],head[1],"pYYS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pzze = cutdata(data[1],head[1],"pZZS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pxye = cutdata(data[1],head[1],"pXYS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pxze = cutdata(data[1],head[1],"pXZS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pyze = cutdata(data[1],head[1],"pYZS0",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
#=
X, Z, Pxxi = cutdata(data[1],head[1],"pXXS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pyyi = cutdata(data[1],head[1],"pYYS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pzzi = cutdata(data[1],head[1],"pZZS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pxyi = cutdata(data[1],head[1],"pXYS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pxzi = cutdata(data[1],head[1],"pXZS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
X, Z, Pyzi = cutdata(data[1],head[1],"pYZS1",cut='y',cutPlaneIndex=cI,
	plotrange=plotrange)
=#

x  = X[:,1]
z  = Z[1,:]

# [A/m^2]
#Jy = @. (qi*ρi*Uyi+qe*ρe*Uye)*1e3/mp*1e6
#Jz = @. (qi*ρi*Uzi+qe*ρe*Uze)*1e3/mp*1e6
const q = 1.6021765e-19 # [C]
# [/cc] --> [/mc], [km] --> [m]
#Jx = @. (qi*ρi/mi*Uxi+qe*ρe/me*Uxe)*1e9*q
#Jy = @. (qi*ρi/mi*Uyi+qe*ρe/me*Uye)*1e9*q
#Jz = @. (qi*ρi/mi*Uzi+qe*ρe/me*Uze)*1e9*q
Jx = @. qi*ρi/mi*Uxi+qe*ρe/me*Uxe
Jy = @. qi*ρi/mi*Uyi+qe*ρe/me*Uye
Jz = @. qi*ρi/mi*Uzi+qe*ρe/me*Uze

# Normalized quantities
fig, ax = plt.subplots(9,2,figsize=(9.0,10.0))
c = Vector{PyObject}(undef,length(ax))
axin = Vector{PyObject}(undef,length(ax))
for i in 1:length(ax)
	axin[i] = inset_axes(ax[i],
      width="5%",  # width = 5% of parent_bbox width
      height="100%",  # height : 50%
      loc="lower left",
      bbox_to_anchor=(1.02, 0., 1.0, 1.0),
      bbox_transform=ax[i].transAxes,
      borderpad=0,)
   axin[i].tick_params(axis="y", direction="in")

   ax[i].tick_params(which="both", direction="in")
end

# Set plotting parameters
levels = 40
plt.set_cmap("seismic")
#vPos, vPos2 = (0.4, 0.8), (0.28,0.8)
vPos, vPos2 = (0.8, 0.8), (0.65,0.8)
lPos = (-0.1, 0.92)
yPos = (-0.17,0.35)

labels = [L"B_z", L"B_y", L"E_x", L"v_{iy}", L"v_{iz}", L"v_{ex}", L"v_{ey}",
	L"v_{ez}", L"\rho_i", L"J_x", L"J_y", L"J_z", L"(E+v_i\times B)_x",
	L"(E+v_e \times B)_x", L"(E+v_i\times B)_y", L"(E+v_e\times B)_y",
	L"D_{ng}", L"D_e"]

vm = ones(18)
const ϵ = 0.05 # allow room for extreme colors
vm[1] = max(abs.(extrema(Bz./B₀))...)+ϵ
vm[2] = max(abs.(extrema(By./B₀))...)+ϵ
vm[3] = max(abs.(extrema(Ex./E₀))...)+ϵ
vm[4] = max(abs.(extrema(Uyi./vAlfven))...)+ϵ
vm[5] = max(abs.(extrema(Uzi./vAlfven))...)+ϵ
vm[6] = max(abs.(extrema(Uxe./vAlfven))...)+ϵ
vm[7] = max(abs.(extrema(Uye./vAlfven))...)+ϵ
vm[8] = max(abs.(extrema(Uze./vAlfven))...)+ϵ

vm[10] = max(abs.(extrema(Jx./J₀))...)+ϵ
vm[11] = max(abs.(extrema(Jy./J₀))...)+ϵ
vm[12] = max(abs.(extrema(Jz./J₀))...)+ϵ
vm[13] = max(abs.(extrema((Ex.+Uyi.*Bz.-Uzi.*By)./E₀))...)+ϵ
vm[14] = max(abs.(extrema((Ex.+Uye.*Bz.-Uze.*By)./E₀))...)+ϵ
vm[15] = max(abs.(extrema((Ey.+Uzi.*Bx.-Uxi.*Bz)./E₀))...)+ϵ
vm[16] = max(abs.(extrema((Ey.+Uze.*Bx.-Uxe.*Bz)./E₀))...)+ϵ
#=
zstart = collect(range(z[10],stop=z[end-10],length=10))
xstart = fill(-1.92,size(zstart))
append!(zstart, collect(range(-0.3,0.4,length=3)))
append!(xstart, fill(-1.82,3))
=#
seeds = select_seeds(x[10:end-10],z[10:end-10]; nSeed=5)
xstart, zstart = seeds[1,:], seeds[2,:]
append!(xstart, [-1.9, -1.9, -1.95, -1.97, -1.95, -1.9, -1.95, -1.95, -1.95])
append!(zstart, [-0.4, -0.5, -0.4, 0.3, 0.4, 0.2, 0.05, -0.1, -0.2])

xl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
zl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
for i = 1:length(xstart)
   xs,zs = xstart[i],zstart[i]
   xl[i], zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds=0.03, maxstep=10000,
   gridType="ndgrid")
end

# Bz
c[1] = ax[1].contourf(Z,X,Bz./B₀,levels, norm=DN(0), vmin=-vm[1], vmax=vm[1])

# By
c[2] = ax[2].contourf(Z,X,By./B₀,levels, norm=DN(0), vmin=-vm[2], vmax=vm[2]) #

# Ex
c[3] = ax[3].contourf(Z,X,Ex./E₀,levels, norm=DN(0), vmin=-vm[3], vmax=vm[3])

# Uyi
c[4] = ax[4].contourf(Z,X,Uyi./vAlfven,levels, norm=DN(0), vmin=-vm[4], vmax=vm[4])

# Uzi
c[5] = ax[5].contourf(Z,X,Uzi./vAlfven,levels, norm=DN(0))

# Uxe
c[6] = ax[6].contourf(Z,X,Uxe./vAlfven,levels, norm=DN(0), vmin=-vm[6], vmax=vm[6])

# Uye
c[7] = ax[7].contourf(Z,X,Uye./vAlfven,levels, norm=DN(0), vmin=-vm[7], vmax=vm[7])

# Uze
c[8] = ax[8].contourf(Z,X,Uze./vAlfven,levels, norm=DN(0), vmin=-vm[8], vmax=vm[8])

# ρi
c[9] = ax[9].contourf(Z,X,ρi./ρ₀,levels, cmap="jet")

# Jx
c[10] = ax[10].contourf(Z,X,Jx./J₀,levels, norm=DN(0), vmin=-vm[10], vmax=vm[10])

# Jy
c[11] = ax[11].contourf(Z,X,Jy./J₀,levels, norm=DN(0), vmin=-vm[11], vmax=vm[11])

# Jz
c[12] = ax[12].contourf(Z,X,Jz./J₀,levels, norm=DN(0), vmin=-vm[12], vmax=vm[12])

# Deviation from ideal MHD
c[13] = ax[13].contourf(Z,X, (Ex.+Uyi.*Bz.-Uzi.*By)./E₀, levels, norm=DN(0),
   vmin=-vm[13], vmax=vm[13])

# Deviation from Hall MHD
c[14] = ax[14].contourf(Z,X, (Ex.+Uye.*Bz.-Uze.*By)./E₀, levels, norm=DN(0),
   vmin=-vm[14], vmax=vm[14])

# Deviation from ideal MHD
c[15] = ax[15].contourf(Z,X, (Ey.+Uzi.*Bx.-Uxi.*Bz)./E₀, levels, norm=DN(0),
   vmin=-vm[15], vmax=vm[15])

# Deviation from Hall MHD
c[16] = ax[16].contourf(Z,X, (Ey.+Uze.*Bx.-Uxe.*Bz)./E₀, levels, norm=DN(0),
   vmin=-vm[16], vmax=vm[16])

# non-gyrotropy index Dng (for electron, not for electron+ion!)
Dng = @. 2*√(Pxye^2 + Pxze^2 + Pyze^2) / (Pxxe + Pyye + Pzze)

c[17] = ax[17].contourf(Z,X,Dng, levels, cmap="jet")

# Dissipation measure De
Dₑ = @. (Jx*(Ex + Uye*Bz - Uze*By) +
		Jy*(Ey + Uze*Bx - Uxe*Bz) +
		Jz*(Ez + Uxe*By - Uye*Bx) -
		(ρi/mi - ρe/me)*(Uxe*Ex + Uye*Ey + Uze*Ez)) / (J₀*B₀*vAlfven)

vm[18] = max(abs.(Dₑ)...)+ϵ

c[18] = ax[18].contourf(Z,X,Dₑ, levels, norm=DN(0), vmin=-vm[18], vmax=vm[18])

for i in 1:length(ax)
   #.ax.locator_params(nbins=5) does not work together with norm(0)!
   cb = colorbar(c[i], cax=axin[i])
   cb.ax.tick_params(labelsize=6)
   cb.ax.locator_params(nbins=5)
   if i in (1,9,17) #(1,4,17)
	   ax[i].annotate(labels[i], xy=vPos, xycoords="axes fraction",color="w")
   elseif i in (13,14,15,16)
      ax[i].annotate(labels[i], xy=vPos2, xycoords="axes fraction")
   else
      ax[i].annotate(labels[i], xy=vPos, xycoords="axes fraction")
   end
   ax[i].annotate("($('a'+i-1))", xy=lPos, xycoords="axes fraction")
   i ≤ length(ax)/2 &&
   ax[i].annotate(L"x [R_G]", xy=yPos, xycoords="axes fraction",rotation=90)
   i % (length(ax)/2) == 0 && ax[i].set_xlabel(L"z [R_G]")
   if i < length(ax)-1
      [ax[i].plot(zl[j],xl[j],"-",color="k",lw=0.4) for j in 1:length(xstart)]
   end
   ax[i].contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
   ax[i].set_aspect("equal", "box")
   ax[i].invert_yaxis()
   i % (length(ax)/2) != 0 && ax[i].axes.xaxis.set_ticklabels([])
   ax[i].tick_params(which="both",top=true, right=true)
   ax[i].minorticks_on()
end

fig.subplots_adjust(wspace=0.02, hspace=0.07)
#tight_layout()


# T∥ and T⟂
# vec0 = (x0, y0, z0) = (0,1,0)
# vec1 = (x1, y1, z1) = (bx, by, bz)/b
# vec2 = (x2, y2, y2) = vec0 x vec1 / norm
# vec3 = (x3, y3, z3) = vec1 x vec2

#=
B = @. √(Bx^2 + By^2 + Bz^2)
x0, y0, z0 = 0.0, 1.0, 0.0
x1, y1, z1 = Bx./B, By./B, Bz./B
x2, y2, z2 = @. y0*z1-z0*y1, z0*x1-x0*z1, x0*y1-y0*x1
z3 = @. sqrt(x2^2 + y2^2 + z2^2)
x2, y2, z2 = @. x2/z3, y2/z3, z2/z3
x3, y3, z3 = @. y1*z2-z1*y2, z1*x2-x1*z2, x1*y2-y1*x2

P11e = @. Pxxe*x1^2+Pyye*y1^2+Pzze*z1^2+2*(Pxye*x1*y1+Pxze*x1*z1+Pyze*y1*z1)
P22e = @. Pxxe*x2^2+Pyye*y2^2+Pzze*z2^2+2*(Pxye*x2*y2+Pxze*x2*z2+Pyze*y2*z2)
P33e = @. Pxxe*x3^2+Pyye*y3^2+Pzze*z3^2+2*(Pxye*x3*y3+Pxze*x3*z3+Pyze*y3*z3)
#P12e = @. Pxxe*x1*x2+Pyye*y1*y2+Pzze*z1*z2+Pxye*(x1*y2+y1*x2)+Pxze*(x1*z2+z1*x2)+Pyze*(y1*z2+z1*y2)
#P13e = @. Pxxe*x1*x3+Pyye*y1*y3+Pzze*z1*z3+Pxye*(x1*y3+y1*x3)+Pxze*(x1*z3+z1*x3)+Pyze*(y1*z3+z1*y3)
#P23e = @. Pxxe*x2*x3+Pyye*y2*y3+Pzze*z2*z3+Pxye*(x2*y3+y2*x3)+Pxze*(x2*z3+z2*x3)+Pyze*(y2*z3+z2*y3)

Tpar = P11e./ρe.*me
Tperp = @. √(P22e^2 + P33e^2)./ρe.*me

figure()
ax1 = subplot(211)
contourf(Z,X, Tpar./T₀, levels, cmap="viridis"); colorbar()
clim([0.0,6.0])
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
ax1.invert_yaxis()
title(L"T_e\parallel")
ax1.axes.xaxis.set_ticklabels([])
ax2 = subplot(212)
contourf(Z,X, Tperp./T₀, levels, cmap="viridis"); colorbar()
clim([0.0,6.0])
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
ax2.invert_yaxis()
title(L"T_e\perp")

#figure()
#contourf(X,Z,Tpar./Tperp, levels); colorbar()
=#

#=
figure()
ax1 = subplot(111)
ax1.invert_yaxis()
contourf(Z,X,Uxe./vAlfven,levels, norm=DN(0))
colorbar()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
=#
