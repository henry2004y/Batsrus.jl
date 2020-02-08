# Plotting PIC variables from box outputs
#
# Hongyang Zhou, hyzhou@umich.edu 02/06/2020

using VisAna, PyPlot

#include("/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia/src/constants.jl")

dir = "/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia"
fnameField = "3d_var_region0_0_t00001640_n00020369.out"

head, data = readdata(fnameField, dir=dir)

#=
pxxe_ = findfirst(x->x=="pXXS0",head[1][:wnames])
pyye_ = findfirst(x->x=="pYYS0",head[1][:wnames])
pzze_ = findfirst(x->x=="pZZS0",head[1][:wnames])
pxye_ = findfirst(x->x=="pXYS0",head[1][:wnames])
pxze_ = findfirst(x->x=="pXZS0",head[1][:wnames])
pyze_ = findfirst(x->x=="pYZS0",head[1][:wnames])
pxxi_ = findfirst(x->x=="pXXS1",head[1][:wnames])
pyyi_ = findfirst(x->x=="pYYS1",head[1][:wnames])
pzzi_ = findfirst(x->x=="pZZS1",head[1][:wnames])
pxyi_ = findfirst(x->x=="pXYS1",head[1][:wnames])
pxzi_ = findfirst(x->x=="pXZS1",head[1][:wnames])
pyzi_ = findfirst(x->x=="pYZS1",head[1][:wnames])
=#

plotrange = [-2.05, -1.75, -0.5, 0.5]
#plotrange=[-Inf, Inf, -Inf, Inf]
cI = 129

# E: [μV/m]
# B: [nT]
# V: [km/s]
# ρ: [amu/cc]
# P: [nPa]

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

x  = X[:,1]
z  = Z[1,:]

# [A/m^2]
#Jy = @. (qi*ρi*Uyi+qe*ρe*Uye)*1e3/mp*1e6
#Jz = @. (qi*ρi*Uzi+qe*ρe*Uze)*1e3/mp*1e6

me = head[1][:eqpar][1]
qe = head[1][:eqpar][2]
mi = head[1][:eqpar][3]
qi = head[1][:eqpar][4]
const vAlfven = 253. # reference Alfven velocity, [km/s]
const B₀ = √((-10.)^2+(-6.)^2+(-86.)^2)
const E₀ = 140.0*B₀ # [μV/m]
const ρ₀ = 56.0     # [amu/cc]

Jx = @. qi*ρi/mi*Uxi+qe*ρe/me*Uxe
Jy = @. qi*ρi/mi*Uyi+qe*ρe/me*Uye
Jz = @. qi*ρi/mi*Uzi+qe*ρe/me*Uze

#=
qi = q
qe = -q
mS0 = 0.14
mS1 = 14
qi = qS1 = 1.0
qe = qS0 =  -1.0
'n0'       , '{rhos0}/mS0'
'n1'       , '{rhos1}/mS1'
'jpz', 'qi*{n1}*{uzs1}+qe*{n0}*{uzs0}'
=#


# Should I normalize all the values?
figure(figsize=(10,8))

# Set plotting parameters
levels = 40
plt.set_cmap("seismic")
labelPos = (0.4, 0.05)
zstart = collect(range(z[10],stop=z[end-10],length=8))
xstart = fill(-2.0,size(zstart))

zstart = append!(zstart, collect(range(-0.21,0.31,length=3)))
xstart = append!(xstart, fill(-1.82,3))

f = subplot(8,2,1)
contourf(Z,X,Bz./B₀,levels) #
colorbar()
#streamplot(Z,X,Bz,Bx, color="k",density=0.5)

xl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
zl = [Vector{Float32}(undef,0) for _ in 1:length(xstart)]
for i = 1:length(xstart)
   xs,zs = xstart[i],zstart[i]
   xl[i], zl[i] = trace2d_rk4(Bx, Bz, xs, zs, x, z, ds=0.02, maxstep=10000, gridType="ndgrid")
   plot(zl[i],xl[i],"-",color="k")
end

plt.gca().invert_yaxis()
axis("scaled")
f.axes.xaxis.set_ticklabels([])
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"B_z", xy=labelPos, xycoords="axes fraction",color="w")


f = subplot(8,2,3)
contourf(Z,X,By./B₀,levels) #
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"B_y", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])


f = subplot(8,2,5)
contourf(Z,X,Ex./E₀,levels)
#colorbar()
#colorbar(boundaries=vplot, ticks=vdisp)
colorbar().ax.tick_params(axis="y", direction="in")
#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
#b = plt.colorbar(f, cax=[0.8, 0.1, 0.03, 0.8])
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"E_x", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])


f = subplot(8,2,7)
contourf(Z,X,Uzi./vAlfven,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{iz}", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,9)
contourf(Z,X,Uyi./vAlfven,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{iy}", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,11)
contourf(Z,X,Uze./vAlfven,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{ez}", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,13)
contourf(Z,X,Uye./vAlfven,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{ey}", xy=labelPos, xycoords="axes fraction",color="w")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,15)
contourf(Z,X,ρi./ρ₀,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"\rho_i", xy=labelPos, xycoords="axes fraction",color="w")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])


f = subplot(8,2,2)
contourf(Z,X,Jz,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"J_z", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,4)
contourf(Z,X,Jy,levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"J_y", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,6)
contourf(Z,X, (Ex.+Uyi.*Bz.-Uzi.*By)./E₀, levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_i\times B)_x", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,8)
contourf(Z,X, (Ex.+Uye.*Bz.-Uze.*By)./E₀, levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_e \times B)_x", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,10)
contourf(Z,X, (Ey.+Uzi.*Bx.-Uxi.*Bz)./E₀, levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_i\times B)_y", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(8,2,12)
contourf(Z,X, (Ey.+Uze.*Bx.-Uxe.*Bz)./E₀, levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_e\times B)_y", xy=labelPos, xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

# non-gyrotropy index Dng
Dng = @. √(Pxye^2 + Pxze^2 + Pyze^2 + Pxyi^2 + Pxzi^2 + Pyzi^2) /
	(Pxxe + Pyye + Pzze + Pxxi + Pyyi + Pzzi)

f = subplot(8,2,14)
contourf(Z,X,Dng, levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"D_{ng}", xy=labelPos, xycoords="axes fraction",color="w")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])


# Dissipation measure De
# The unit is not right!!!
Dₑ = @. Jx*(Ex + Uye*Bz - Uze*By) +
		Jy*(Ey + Uze*Bx - Uxe*Bz) +
		Jz*(Ez + Uxe*By - Uye*Bx) -
		(ρi/mi - ρe/me)*(Uxe*Ex + Uye*Ey + Uze*Ez)

f = subplot(8,2,16)
contourf(Z,X,Dₑ, levels)
colorbar().ax.tick_params(axis="y", direction="in")
[plot(zl[i],xl[i],"-",color="k") for i in 1:length(xstart)]
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"D_e", xy=labelPos, xycoords="axes fraction")
#f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

#tight_layout()
