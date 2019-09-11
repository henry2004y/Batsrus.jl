using PyCall, PyPlot

## Example 1
py"""
import numpy as np
def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
"""

np = pyimport("numpy")
scipy = pyimport("scipy")
interpolate = pyimport("scipy.interpolate")
griddata = interpolate.griddata

s = pybuiltin(:slice)
XY = get(np.mgrid, (s(0,1,100im), s(0,1,200im)))
grid_x, grid_y = XY[1,:,:], XY[2,:,:]

points = np.random.rand(1000, 2)
values = py"func"(points[:,1], points[:,2])

grid_z0 = griddata(points, values, (grid_x, grid_y), method="nearest")
grid_z1 = griddata(points, values, (grid_x, grid_y), method="linear")
grid_z2 = griddata(points, values, (grid_x, grid_y), method="cubic")

plt.subplot(221)
plt.imshow(py"func"(grid_x, grid_y), extent=(0,1,0,1), origin="lower")
plt.plot(points[:,1], points[:,2], "k.", ms=1)
plt.title("Original")
plt.subplot(222)
plt.imshow(grid_z0, extent=(0,1,0,1), origin="lower")
plt.title("Nearest")
plt.subplot(223)
plt.imshow(grid_z1, extent=(0,1,0,1), origin="lower")
plt.title("Linear")
plt.subplot(224)
plt.imshow(grid_z2, extent=(0,1,0,1), origin="lower")
plt.title("Cubic")
plt.gcf().set_size_inches(6, 6)
plt.show()



## Example 2
py"""
def func(x, y):
    return x + y
"""

s = pybuiltin(:slice)
XY = get(np.mgrid, (s(0,1,5im),s(0,1,5im)))
grid_x, grid_y = XY[1,:,:], XY[2,:,:]
points = np.random.rand(25, 2)
values = py"func"(points[:, 1], points[:, 2])

griddata(points, values, (grid_x, grid_y), method="nearest")
griddata(points, values, (grid_x, grid_y), method="cubic")

using PyPlot
plt.plot(points[:,1], points[:,2], "k.", ms=1)

## test 1
py"""
def func(x, y, z):
    return x + y + z
"""
np = pyimport("numpy")
scipy = pyimport("scipy")
interpolate = pyimport("scipy.interpolate")
griddata = interpolate.griddata

s = pybuiltin(:slice)
XYZ = get(np.mgrid, (s(0,1,10im), s(0,1,20im), s(0,1,10im)))
grid_x, grid_y, grid_z = XYZ[1,:,:,:], XYZ[2,:,:,:], XYZ[3,:,:,:]

points = np.random.rand(1000, 3)
values = py"func"(points[:,1], points[:,2], points[:,3])

grid_z0 = griddata(points, values, (grid_x, grid_y, grid_z), method="nearest")
grid_z1 = griddata(points, values, (grid_x, grid_y, grid_z), method="linear")
grid_z2 = griddata(points, values, (grid_x, grid_y, grid_z), method="cubic")
