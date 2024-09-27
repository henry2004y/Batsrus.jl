# Examples

## IDL format output loader

- Read data

```julia
file = "1d_bin.out";
bd = load(file);
bd = load(file, verbose=true);
bd = load(file, npict=1);
```

- 3D structured spherical coordinates

```julia
file = "3d_structured.out";
bd = load(file, verbose=false);
```

- log file

```julia
logfilename = "shocktube.log";
head, data = readlogdata(logfilename)
```

### Data Extraction

- Checking variable range

```julia
get_var_range(bd, "rho")
```

- Raw variables

```julia
ρ = getvar(bd, "rho")
bd["rho"]
```

- Derived variables

```julia
v = getvars(bd, ["Bx", "By", "Bz"])
Bmag = bd["Bmag"]
```

- Extracting data at a given location

```julia
loc = Float32[0.0, 0.0] # The type determines the output type
d = interp1d(bd, "rho", loc)
```

- Extracting data along a given line

```julia
point1 = Float32[-10.0, -1.0]
point2 = Float32[10.0, 1.0]
w = interp1d(bd, "rho", point1, point2)
```

Here is a full list of predefined derived quantities:

| Derived variable name | Meaning                          | Required variable |
|-----------------------|----------------------------------|-------------------|
| B2                    | magnetic field magnitude squared | Bx, By, Bz        |
| E2                    | electric field magnitude squared | Ex, Ey, Ez        |
| U2                    | velocity magnitude squared       | Ux, Uy, Uz        |
| Bmag                  | magnetic field magnitude         | Bx, By, Bz        |
| Emag                  | electric field magnitude         | Ex, Ey, Ez        |
| Umag                  | velocity magnitude               | Ux, Uy, Uz        |
| B                     | magnetic field vector            | Bx, By, Bz        |
| E                     | electric field vector            | Ex, Ey, Ez        |
| U                     | velocity vector                  | Ux, Uy, Uz        |

## Output format conversion

We can convert 2D/3D BATSRUS outputs `*.dat` to VTK formats. It uses the VTK XML format writer [writeVTK](https://github.com/jipolanco/WriteVTK.jl) to generate files for Paraview and Tecplot. The default converted filename is `out.vtu`.

- ASCII Tecplot file (supports both `tec` and `tcp`) and binary Tecplot file (set `DOSAVETECBINARY=TRUE` in BATSRUS `PARAM.in`):

```julia
file = "x=0_mhd_1_n00000050.dat"
convertTECtoVTU(file)
```

- 3D structured IDL file (`gridType=1` returns rectilinear `vtr` file, `gridType=2` returns structured `vts` file):

```julia
file = "3d_structured.out"
convertIDLtoVTK(file, gridType=1)
```

- 3D unstructured IDL file together with header and tree file:

```julia
filetag = "3d_var_1_n00002500"
convertIDLtoVTK(filetag)
```

!!! note
    The file suffix should not be provided for this to work correctly!

- Multiple files:

```julia
dir = "./"
filenames = filter(file -> startswith(file, "3d") && endswith(file, ".dat"), readdir(dir))
filenames = dir .* filenames

for filename in filenames
   convertTECtoVTU(filename, filename[1:end-4])
end
```

* Processing multiple files with threads in parallel:

```julia
dir = "./"
filenames = filter(file -> startswith(file, "3d") && endswith(file, ".dat"), readdir(dir))
filenames = dir .* filenames

Threads.@threads for filename in filenames
   println("filename=$filename")
   convertTECtoVTU(filename, filename[1:end-4])
end
```

More examples can be found in [examples](https://github.com/henry2004y/Batsrus.jl/tree/master/examples).

## HDF format output loader

```julia
filename = "3d__var_1_n00006288.h5"
file = BatsrusHDF5Uniform(filename)
```

### Field extraction

Variable `var` can be extracted in the whole domain:

```julia
var, (xl_new, yl_new, zl_new), (xu_new, yu_new, zu_new) = extract_var(file, "bx")
```

where `(xl_new, yl_new, zl_new)` and `(xu_new, yu_new, zu_new)` return the lower and upper bound, respectively.

Variables within a box region can be extracted as following:

```julia
var, (xl_new, yl_new, zl_new), (xu_new, yu_new, zu_new) =
   extract_var(file, "bx"; xmin, xmax, ymin, ymax, zmin, zmax)
```

## Data visualization

We provide plot recipes for Plots.jl, Makie.jl, and wrappers for PyPlot.jl.

The recipes for Plots.jl and Makie.jl will work on all kinds of plots given the correct dimensions, e.g.

```julia
using Plots

plot(bd, "p")
contourf(bd, "Mx", xlabel="x")
```

See the official documentation for Plots.jl for more information.

On the other hand, most common 1D and 2D plotting functions are wrapped over their Matplotlib equivalences through PyPlot.jl.
To trigger the wrapper, `using PyPlot`.
Check out the documentation for more details.

### Quick exploration of data

A general `plotdata` function is provided for quick visualizations using Matplotlib.

- 1D binary

```julia
plotdata(bd, "p", plotmode="line")
plotdata(bd, "p", plotmode="linegrid")
```

- 2D Cartesian (structured)

```julia
plotdata(bd, "p bx;by", plotmode="contbar streamover")
plotdata(bd, "p bx;by", plotmode="contbar quiverover")
plotdata(bd, "p bx;by", plotmode="contbar streamover", density=2.0)
plotdata(bd, "p", plotmode="grid")
plotdata(bd, "p", plotmode="contbar", plotrange=[-50., 50., -1., 1.])
plotdata(bd, "p", plotmode="contbar")
plotdata(bd, "p", plotmode="contbarlog")
plotdata(bd, "p", plotmode="surfbar")
```

- 2D unstructured

```julia
plotdata(bd, "rho", plotmode="contbar")
plotdata(bd, "rho", plotmode="trimesh")
plotdata(bd, "rho", plotmode="tricont")
```

- 2D structured spherical coordinates

```julia
plotdata(bd, "rho", plotmode="contbar")
```

- 3D box

```julia
plotdata(bd, "bx", plotmode="contbar", dir="y", sequence=1, level=20)
plotdata(bd, "bx", plotmode="contbar", dir="y", plotrange=[-1.4,-1.1,0.70,0.78])
plt.axis("scaled")

subplot(2,2,(1,3))
cutplot(bd, "Ex"; dir="y", sequence=128, plotrange)
```

#### Finding indexes

To get the index of a certain quantity, e.g. electron number density

```julia
ρe_= findfirst(x->x=="rhoS0", bd.head.wnames)
```

### Multiple dispatch for Matplotlib functions

Using the same plotting functions as in Matplotlib is allowed, and actually recommended.
Some plotting functions can be directly called as shown below, which allows for more control from the user.
`using PyPlot` to import the full capability of the package, etc. adding colorbar, changing line colors, setting colorbar range with `clim`.

- line plot

```julia
plot(bd, "p", linewidth=2, color="green")
c = plot(bd, "p")
plt.setp(c, linestyle="--", linewidth=2);
```

- scatter plot

```julia
scatter(bd, "p")
```

- contour

```julia
contour(bd, "p")
```

- filled contour

```julia
contourf(bd, "p")
contourf(bd, "p", levels, plotrange=[-10,10,-Inf,Inf], plotinterval=0.1)
```

- surface plot

```julia
plot_surface(bd, "p")
```

- triangle surface plot

```julia
plot_trisurf(bd, "p")
```

- triangle filled contour plot

```julia
tricontourf(bd, "p")
```

- streamline

```julia
streamplot(bd, "bx;bz")
streamplot(bd, "bx;bz", density=2.0, color="k", plotinterval=1.0, plotrange=[-10,10,-Inf,Inf])
```

- quiver (currently only for Cartesian grid)

```julia
quiver(bd, "ux;uy", stride=50)
```

- streamline + contourf

```julia
file = "y.out"
bd = load(file)

DN = matplotlib.colors.DivergingNorm
cmap = matplotlib.cm.RdBu_r

contourf(bd, "uxS0", 50; plotrange=[-3,3,-3,3], plotinterval=0.05, norm=DN(0), cmap)
colorbar()
streamplot(bd, "uxS0;uzS0"; density=2.0, color="g", plotrange=[-3,3,-3,3])
xlabel("x"); ylabel("y"); title("Ux [km/s]")

contourf(bd, "uxS0", 50; plotinterval=0.05, norm=DN(0), cmap)
colorbar()
axis("scaled")
xlabel("x"); ylabel("y"); title("uxS0")
```

### Tracing

The built-in `streamplot` function in Matplotlib is not satisfactory for accurately tracing. Instead we recommend [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) for tracing fieldlines and streamlines.

An example of tracing in a 2D cut and plot the field lines over contour:

```julia
file = "test/y=0_var_1_t00000000_n00000000.out"
bd = load(file)

bx = bd.w[:,:,5]
bz = bd.w[:,:,7]
x  = bd.x[:,1,1]
z  = bd.x[1,:,2]

seeds = select_seeds(x, z; nSeed=100) # randomly select the seeding points

for i in 1:size(seeds)[2]
   xs = seeds[1,i]
   zs = seeds[2,i]
   # Tracing in both direction. Check the document for more options.
   x1, z1 = trace2d_eul(bx, bz, xs, zs, x, z, ds=0.1, maxstep=1000, gridType="ndgrid")
   plot(x1, z1, "--")
end
axis("equal")
```

Currently the `select_seeds` function uses pseudo random number generator that produces the same seeds every time.
