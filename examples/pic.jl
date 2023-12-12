# Plotting PIC variables from box outputs in normalized units.
#
# Hongyang Zhou, hyzhou@umich.edu 02/06/2020

using Batsrus, PyPlot

# For precise colorbar control
using PyCall
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
inset_axes = inset_locator.inset_axes
# For specifying the zero point in the colorbar
DN = matplotlib.colors.DivergingNorm

dir = "."
file = "3d_var_region0_0_t00001640_n00020369.out"

data = load(joinpath(dir, file))

me = data.head.eqpar[1]
qe = data.head.eqpar[2]
mi = data.head.eqpar[3]
qi = data.head.eqpar[4]
const kB = 1.38064852e-23 # [m^2 kg s^-2 K^-1]
const q = 1.6021765e-19 # [C]
const vAlfven = 253. # reference Alfven velocity, [km/s]
const B₀ = √((-10.)^2+(-6.)^2+(-86.)^2)
const E₀ = vAlfven*B₀ # [μV/m]
const ρ₀ = 56.0     # [amu/cc]
const J₀ = 4.0*vAlfven # Actual normalization unit: q*4.0*vAlfven
#const T₀ = vAlfven^2
const T₀ = 0.2/4 # Pe/n₀

plotrange = [-2.05, -1.75, -0.5, 0.5]
#plotrange=[-Inf, Inf, -Inf, Inf]
cI = 129 # plane cut index

##
# ρ: [amu/cc]
# B: [nT]
# E: [μV/m]
# V: [km/s]
# P: [nPa]

X, Z, ρe = cutdata(data, "rhoS0", dir="y", sequence=cI, plotrange)
X, Z, ρi = cutdata(data, "rhoS1", dir="y", sequence=cI, plotrange)
X, Z, Bx = cutdata(data, "Bx",    dir="y", sequence=cI, plotrange)
X, Z, By = cutdata(data, "By",    dir="y", sequence=cI, plotrange)
X, Z, Bz = cutdata(data, "Bz",    dir="y", sequence=cI, plotrange)
X, Z, Ex = cutdata(data, "Ex",    dir="y", sequence=cI, plotrange)
X, Z, Ey = cutdata(data, "Ey",    dir="y", sequence=cI, plotrange)
X, Z, Ez = cutdata(data, "Ez",    dir="y", sequence=cI, plotrange)
X, Z, Uxe= cutdata(data, "uxS0",  dir="y", sequence=cI, plotrange)
X, Z, Uye= cutdata(data, "uyS0",  dir="y", sequence=cI, plotrange)
X, Z, Uze= cutdata(data, "uzS0",  dir="y", sequence=cI, plotrange)
X, Z, Uxi= cutdata(data, "uxS1",  dir="y", sequence=cI, plotrange)
X, Z, Uyi= cutdata(data, "uyS1",  dir="y", sequence=cI, plotrange)
X, Z, Uzi= cutdata(data, "uzS1",  dir="y", sequence=cI, plotrange)

X, Z, Pxxe = cutdata(data, "pXXS0", dir="y", sequence=cI, plotrange)
X, Z, Pyye = cutdata(data, "pYYS0", dir="y", sequence=cI, plotrange)
X, Z, Pzze = cutdata(data, "pZZS0", dir="y", sequence=cI, plotrange)
X, Z, Pxye = cutdata(data, "pXYS0", dir="y", sequence=cI, plotrange)
X, Z, Pxze = cutdata(data, "pXZS0", dir="y", sequence=cI, plotrange)
X, Z, Pyze = cutdata(data, "pYZS0", dir="y", sequence=cI, plotrange)
#=
X, Z, Pxxi = cutdata(data, "pXXS1", dir="y", sequence=cI, plotrange)
X, Z, Pyyi = cutdata(data, "pYYS1", dir="y", sequence=cI, plotrange)
X, Z, Pzzi = cutdata(data, "pZZS1", dir="y", sequence=cI, plotrange)
X, Z, Pxyi = cutdata(data, "pXYS1", dir="y", sequence=cI, plotrange)
X, Z, Pxzi = cutdata(data, "pXZS1", dir="y", sequence=cI, plotrange)
X, Z, Pyzi = cutdata(data, "pYZS1", dir="y", sequence=cI, plotrange)
=#

x  = X[:,1]
z  = Z[1,:]

# [A/m^2]
#Jy = @. (qi*ρi*Uyi+qe*ρe*Uye)*1e3/mp*1e6
#Jz = @. (qi*ρi*Uzi+qe*ρe*Uze)*1e3/mp*1e6

# [/cc] --> [/mc], [km] --> [m]
#Jx = @. (qi*ρi/mi*Uxi+qe*ρe/me*Uxe)*1e9*q
#Jy = @. (qi*ρi/mi*Uyi+qe*ρe/me*Uye)*1e9*q
#Jz = @. (qi*ρi/mi*Uzi+qe*ρe/me*Uze)*1e9*q
Jx = @. qi*ρi/mi*Uxi+qe*ρe/me*Uxe
Jy = @. qi*ρi/mi*Uyi+qe*ρe/me*Uye
Jz = @. qi*ρi/mi*Uzi+qe*ρe/me*Uze

# Normalized quantities
fig, ax = plt.subplots(10,2,figsize=(8.00,9.00))
c = Vector{PyObject}(undef,length(ax))
axin = Vector{PyObject}(undef,length(ax))
for i in eachindex(ax)
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
vPos, vPos2 = (0.4, 0.76), (0.28,0.78)
#vPos, vPos2 = (0.80, 0.75), (0.62,0.77)
lPos = (-0.12, 0.94)
yPos = (-0.24, 0.33)

labels = [L"B_z", L"B_y", L"E_x", L"v_{iy}", L"v_{iz}", L"v_{ex}", L"v_{ey}",
	L"v_{ez}", L"\rho_i", L"J_x", L"J_y", L"J_z", L"(E+v_i\times B)_x",
	L"(E+v_e \times B)_x", L"(E+v_i\times B)_y", L"(E+v_e\times B)_y",
	L"A\o", L"D_{ng}", L"\sqrt{Q}", L"D_e"]

vm = ones(20)
const ϵ = 0.05 # allow room for extreme colors
vm[1] = max(abs.(extrema(Bz./B₀))...) + ϵ
vm[2] = max(abs.(extrema(By./B₀))...) + ϵ
vm[3] = max(abs.(extrema(Ex./E₀))...) + ϵ
vm[4] = max(abs.(extrema(Uyi./vAlfven))...) + ϵ
vm[5] = max(abs.(extrema(Uzi./vAlfven))...) + ϵ
vm[6] = max(abs.(extrema(Uxe./vAlfven))...) + ϵ
vm[7] = max(abs.(extrema(Uye./vAlfven))...) + ϵ
vm[8] = max(abs.(extrema(Uze./vAlfven))...) + ϵ

vm[10] = max(abs.(extrema(Jx./J₀))...) + ϵ
vm[11] = max(abs.(extrema(Jy./J₀))...) + ϵ
vm[12] = max(abs.(extrema(Jz./J₀))...) + ϵ
vm[13] = max(abs.(extrema((Ex.+Uyi.*Bz.-Uzi.*By)./E₀))...) + ϵ
vm[14] = max(abs.(extrema((Ex.+Uye.*Bz.-Uze.*By)./E₀))...) + ϵ
vm[15] = max(abs.(extrema((Ey.+Uzi.*Bx.-Uxi.*Bz)./E₀))...) + ϵ
vm[16] = max(abs.(extrema((Ey.+Uze.*Bx.-Uxe.*Bz)./E₀))...) + ϵ

seeds = select_seeds(x[10:end-10],z[10:end-10]; nSeed=5)
xstart, zstart = seeds[1,:], seeds[2,:]
append!(xstart, [-1.9, -1.9, -1.95, -1.97, -1.95, -1.9, -1.95, -1.95, -1.95])
append!(zstart, [-0.4, -0.5, -0.4, 0.3, 0.4, 0.2, 0.05, -0.1, -0.2])

xl = [Vector{Float32}(undef,0) for _ in eachindex(xstart)]
zl = [Vector{Float32}(undef,0) for _ in eachindex(xstart)]
for i in eachindex(xstart)
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
#c[6] = ax[6].contourf(Z,X,Uxe./vAlfven,levels, norm=DN(0), vmin=-vm[6], vmax=vm[6])
c[6] = ax[6].contourf(Z,X,Uxe./vAlfven,levels, norm=DN(0, vmin=-vm[6], vmax=vm[6]))

# Uye
c[7] = ax[7].contourf(Z,X,Uye./vAlfven,levels, norm=DN(0), vmin=-vm[7], vmax=vm[7])

# Uze
c[8] = ax[8].contourf(Z,X,Uze./vAlfven,levels, norm=DN(0), vmin=-vm[8], vmax=vm[8])

# ρi
c[9] = ax[9].contourf(Z,X,ρi./ρ₀,levels, cmap="inferno")

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

# agyrotropy measure A 8[Scudder 2008]
B² = @. Bx*Bx + By*By + Bz*Bz
Nxx = @. (By*By*Pzze - 2*By*Bz*Pyze + Bz*Bz*Pyye )/B²
Nxy = @. (-By*Bx*Pzze + By*Bz*Pxze + Bz*Bx*Pyze - Bz*Bz*Pxye)/B²
Nxz = @. (By*Bx*Pyze - By*By*Pxze - Bz*Bx*Pyye + Bz*By*Pxye)/B²
Nyy = @. (Bx*Bx*Pzze - 2*Bx*Bz*Pxze + Bz*Bz*Pxxe)/B²
Nyz = @. (-Bx*Bx*Pyze + Bx*By*Pxze + Bz*Bx*Pxye - Bz*By*Pxxe)/B²
Nzz = @. (Bx*Bx*Pyye - 2*Bx*By*Pxye + By*By*Pxxe)/B²
α = @. Nxx + Nyy + Nzz
β = @. -(Nxy*Nxy + Nxz*Nxz + Nyz*Nyz - Nxx*Nyy - Nxx*Nzz - Nyy*Nzz)
A = @. 2*√(α*α - 4β)/α

c[17] = ax[17].contourf(Z,X, A, levels, cmap="inferno")

# non-gyrotropy measure Dng (for electron, not for electron+ion!) [Aunai 2013]
Dng = @. 2*√(Pxye*Pxye + Pxze*Pxze + Pyze*Pyze) / (Pxxe + Pyye + Pzze)
c[18] = ax[18].contourf(Z,X, Dng, levels, cmap="inferno")

# non-gyrotropy measure Q [Swisdak 2016]
I₁ = @. Pxxe + Pyye + Pzze
I₂ = @. Pxxe*Pyye + Pxxe*Pzze + Pyye*Pzze - Pxye*Pxye - Pyze*Pyze - Pxze*Pxze

Ppar = @. (Bx*Bx*Pxxe + By*By*Pyye + Bz*Bz*Pzze +
	2*(Bx*By*Pxye + Bx*Bz*Pxze + By*Bz*Pyze))/B²
Qsqr = @. √(1 - 4I₂/((I₁ - Ppar)*(I₁ + 3Ppar)))

c[19] = ax[19].contourf(Z,X, Qsqr, levels, cmap="inferno")

# Dissipation measure De
Dₑ = @. (Jx*(Ex + Uye*Bz - Uze*By) +
		Jy*(Ey + Uze*Bx - Uxe*Bz) +
		Jz*(Ez + Uxe*By - Uye*Bx) -
		(ρi/mi - ρe/me)*(Uxe*Ex + Uye*Ey + Uze*Ez)) / (J₀*B₀*vAlfven)

vm[20] = max(abs.(Dₑ)...) + ϵ

c[20] = ax[20].contourf(Z, X, Dₑ, levels, norm=DN(0), vmin=-vm[20], vmax=vm[20])

for i in eachindex(ax)
   #.ax.locator_params(nbins=5) does not work together with norm(0)!
   cb = colorbar(c[i], cax=axin[i])
   cb.ax.tick_params(labelsize=5)
   #cb.ax.locator_params(nbins=5)
   if i in (1,9,17,18,19) #(1,9,17,18,19)
      ax[i].annotate(labels[i], xy=vPos, xycoords="axes fraction", color="w",
         weight="bold")
   elseif i in (13,14,15,16)
      ax[i].annotate(labels[i], xy=vPos2, xycoords="axes fraction")
   else
      ax[i].annotate(labels[i], xy=vPos, xycoords="axes fraction")
   end
   ax[i].annotate("($('a'+i-1))", xy=lPos, xycoords="axes fraction")
   i ≤ length(ax)/2 &&
   ax[i].annotate(L"x [R_G]", xy=yPos, xycoords="axes fraction",rotation=90)
   i % (length(ax)/2) == 0 && ax[i].set_xlabel(L"z [R_G]")
   if i < length(ax)-3
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
contourf(Z,X, Tpar./T₀, levels, cmap="inferno"); colorbar()
clim([0.0,6.0])
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
ax1.invert_yaxis()
title(L"T_e\parallel")
ax1.axes.xaxis.set_ticklabels([])
ax2 = subplot(212)
contourf(Z,X, Tperp./T₀, levels, cmap="inferno"); colorbar()
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
