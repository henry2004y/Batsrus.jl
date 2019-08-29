using PyCall, PyPlot

py"""
import numpy as np
def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
"""

np = pyimport("numpy")
scipy = pyimport("scipy")
griddata = scipy.interpolate.griddata

s = pybuiltin(:slice)
XY = get(np.mgrid, (s(0,100,1), s(0,200,1)))
grid_x, grid_y = XY[1,:,:], XY[2,:,:]

points = np.random.rand(1000, 2)
values = py"func"(points[:,1], points[:,2])

from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method="nearest")
grid_z1 = griddata(points, values, (grid_x, grid_y), method="linear")
grid_z2 = griddata(points, values, (grid_x, grid_y), method="cubic")

plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin="lower")
plt.plot(points[:,0], points[:,1], "k.", ms=1)
plt.title("Original")
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin="lower")
plt.title("Nearest")
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin="lower")
plt.title("Linear")
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin="lower")
plt.title("Cubic")
plt.gcf().set_size_inches(6, 6)
plt.show()