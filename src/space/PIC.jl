# Plotting PIC variables from box outputs
#
# Hongyang Zhou, hyzhou@umich.edu 02/06/2020

using VisAna, PyPlot

include("/Users/hyzhou/Ganymede/scripts/VisAnaJulia/src/constants.jl")

dir = "/Users/hyzhou"
fnameField = "3d_var_region0_0_t00001640_n00020369.out"

head, data = readdata(fnameField, dir=dir)

#=
# B, E, V, J, n, P
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


qi = q
qe = -q

plotrange = [-2.05, -1.75, -0.5, 0.5]

X, Z, ρe = cutdata(data[1],head[1],"rhoS0",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, ρi = cutdata(data[1],head[1],"rhoS1",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Bx = cutdata(data[1],head[1],"Bx",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, By = cutdata(data[1],head[1],"By",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Bz = cutdata(data[1],head[1],"Bz",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Ex = cutdata(data[1],head[1],"Ex",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Ey = cutdata(data[1],head[1],"Ey",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Ez = cutdata(data[1],head[1],"Ez",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Uxe= cutdata(data[1],head[1],"uxS0",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Uye= cutdata(data[1],head[1],"uyS0",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Uze= cutdata(data[1],head[1],"uzS0",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Uxi= cutdata(data[1],head[1],"uxS1",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Uyi= cutdata(data[1],head[1],"uyS1",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)
X, Z, Uzi= cutdata(data[1],head[1],"uzS1",cut='y',cutPlaneIndex=128,
	plotrange=plotrange)

#jx = @. qi*ρi*Uxi+qe*ρe*Uxe
Jy = @. qi*ρi*mp*Uyi+qe*ρe*me*Uye
Jz = @. qi*ρi*mp*Uzi+qe*ρe*me*Uze


figure(figsize=(10,8))

levels = 40
plt.set_cmap("seismic")

f = subplot(7,2,1)
contourf(Z,X,Bz,levels)
plt.gca().invert_yaxis()
axis("scaled")
f.axes.xaxis.set_ticklabels([])
#colorbar(orientation="horizontal")
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"B_z", xy=(0.6, 0.05), xycoords="axes fraction")


f = subplot(7,2,3)
contourf(Z,X,By,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"B_y", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

f = subplot(7,2,5)
contourf(Z,X,Ex,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"E_x", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

f = subplot(7,2,7)
contourf(Z,X,Uzi,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{iz}", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

f = subplot(7,2,9)
contourf(Z,X,Uyi,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{iy}", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

f = subplot(7,2,11)
contourf(Z,X,Uze,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{ez}", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

f = subplot(7,2,13)
contourf(Z,X,Uye,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"v_{ey}", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

f = subplot(7,2,2)
contourf(Z,X,ρi,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"\rho_i", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])
#colorbar(orientation="horizontal")

# JL is so different from Shay et. al!!!
f = subplot(7,2,4)
contourf(Z,X,Jz,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"J_z", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(7,2,6)
contourf(Z,X,Jy,levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"J_y", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(7,2,8)
contourf(Z,X, Ex.+Uyi.*Bz.-Uzi.*By, levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_i\times B)_x", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(7,2,10)
contourf(Z,X, Ex.+Uye.*Bz.-Uze.*By, levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_e \times B)_x", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(7,2,12)
contourf(Z,X, Ey.+Uzi.*Bx.-Uxi.*Bz, levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_i\times B)_y", xy=(0.6, 0.05), xycoords="axes fraction")
f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

f = subplot(7,2,14)
contourf(Z,X, Ey.+Uze.*Bx.-Uxe.*Bz, levels)
axis("scaled")
plt.gca().invert_yaxis()
contour(Z,X,Bz,[0.],colors="k",linestyles="dotted",linewidths=1.)
annotate(L"(E+v_e\times B)_y", xy=(0.6, 0.05), xycoords="axes fraction")
#f.axes.xaxis.set_ticklabels([])
f.axes.yaxis.set_ticklabels([])

#tight_layout()

# non-gyrotropy index Dng


# Dissipation measure De
