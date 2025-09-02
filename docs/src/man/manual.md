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

```julia
bd["rho"][X=-10 .. 10, Y=Near(0.0)]
```

- Extracting data using DimensionalData

We can also use [Selectors](https://rafaqz.github.io/DimensionalData.jl/stable/selectors) from DimensionalData for extracting data. Note that the Selectors need to be imported from Batsrus.jl; alternatively you can simply `using DimensionalData`.

```julia
bd["rho"][X=At(0.0), Y=At(0.0)]
bd["rho"][X=-10 .. 10, Y=Near(0.5)]
bd["rho"][X=-10 .. 10, Y=-0.5 .. 0.5]
```

- Derived variables

We provide utility methods `get_magnitude`, `get_magnitude2`, and `fill_vector_from_scalars` for vector processing:

```julia
Bmag = get_magnitude(bd, :B)
B2 = get_magnitude2(bd, :B)
Bvec = Batsrus.fill_vector_from_scalars(bd, :B)
paniso0 = get_anisotropy(bd, 0)
```

These are built upon `get_vectors`. `fill_vector_from_scalars` is slower than `get_vectors` since it involves additional array allocations. Here is a full list of predefined derived quantities in `get_vectors`:

| Derived variable name | Meaning                          | Required variable |
|-----------------------|----------------------------------|-------------------|
| :B                    | Magnetic field vector            | Bx, By, Bz        |
| :E                    | Electric field vector            | Ex, Ey, Ez        |
| :U                    | Velocity vector                  | Ux, Uy, Uz        |
| :U0                   | Electron velocity vector         | UxS0, UyS0, UzS0  |
| :U1                   | Proton velocity vector           | UxS1, UyS1, UzS1  |

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

Using the same plotting functions as in Matplotlib is allowed, and actually recommended. This takes advantage of multiple dispatch mechanism in Julia.
Some plotting functions can be directly called as shown below, which allows for more control from the user.
`using PyPlot` to import the full capability of the package, etc. adding colorbar, changing line colors, setting colorbar range with `clim`.

For 1D outputs, we can use `plot` or `scatter`.

- line plot

```julia
plot(bd, "p", linewidth=2, color="tab:red", linestyle="--", linewidth=2)
```

- scatter plot

```julia
scatter(bd, "p")
```

For 2D outputs, we can select the following functions:

- `contour`
- `contourf`
- `imshow`
- `pcolormesh`
- `plot_surface`
- `plot_tricontour`
- `plot_tricontourf`
- `plot_trisurf`
- `tripcolor`

with either `quiver` or `streamplot`. By default the linear colorscale is applied. If you want to switch to logarithmic, set argument `colorscale=:log`.

- contour

```julia
contour(bd, "p")
```

- filled contour

```julia
contourf(bd, "p")
contourf(bd, "p"; levels, plotrange=[-10,10,-Inf,Inf], plotinterval=0.1)
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
streamplot(bd, "bx;bz"; density=2.0, color="k", plotinterval=1.0, plotrange=[-10,10,-Inf,Inf])
```

- quiver (currently only for Cartesian grid)

```julia
quiver(bd, "ux;uy"; stride=50)
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

For 3D outputs, we may use `cutplot` for visualizing on a sliced plane, or `streamslice` to plot streamlines on a given slice.

#### Finding indexes

To get the index of a certain quantity, e.g. electron number density

```julia
ρe_= findfirst(x->x=="rhoS0", bd.head.wname)
```

#### Get variable range

```julia
wmin, wmax = get_var_range(bd, var)
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
