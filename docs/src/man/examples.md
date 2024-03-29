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

## Data Extraction

- Raw variables

```julia
ρ = getvar(bd, "rho")
bd["rho"]
```

- Derived variables

```julia
v = getvars(bd, ["Bx", "By", "Bz"])
B = @. sqrt(v["Bx"]^2 + v["By"]^2 + v["Bz"]^2)
```

## Output format conversion

We can convert 2D/3D BATSRUS outputs `*.dat` to VTK formats. It uses the VTK XML format writer [writeVTK](https://github.com/jipolanco/WriteVTK.jl) to generate files for Paraview and Tecplot. The default converted filename is `out.vtu`.

ASCII Tecplot file (supports both `tec` and `tcp`) and binary Tecplot file (set `DOSAVETECBINARY=TRUE` in BATSRUS `PARAM.in`):

```julia
file = "x=0_mhd_1_n00000050.dat"
#file = "3d_ascii.dat"
#file = "3d_bin.dat"
head, data, connectivity = readtecdata(file)
convertTECtoVTU(head, data, connectivity)
```

3D structured IDL file (`gridType=1` returns rectilinear `vtr` file, `gridType=2` returns structured `vts` file):

```julia
file = "3d_structured.out"
convertIDLtoVTK(file, gridType=1)
```

3D unstructured IDL file together with header and tree file:

```julia
filetag = "3d_var_1_n00002500"
convertIDLtoVTK(filetag)
```

!!! note
    The file suffix should not be provided for this to work correctly!

Multiple files:

```julia
using Batsrus, Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef, 0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
tec = readtecdata.(filenames) # head, data, connectivity
for (i, outname) in enumerate(filenames)
   convertTECtoVTU(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])
end
```

If each individual file size is large, consider using:

```julia
using Batsrus, Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef, 0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
for (i, outname) in enumerate(filenames)
   head, data, connectivity = readtecdata(outname)
   convertTECtoVTU(head, data, connectivity, outname[1:end-4])
end
```

Multiple files in parallel:

```julia
using Distributed
@everywhere using Batsrus, Glob

filenamesIn = "cut*.dat"
dir = "."
filenames = Vector{String}(undef, 0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname)
   convertTECtoVTU(head, data, connectivity, outname[1:end-4])
end
```

More examples can be found in [examples](https://github.com/henry2004y/Batsrus.jl/tree/master/examples).

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
using PyPlot
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
# 2D contour
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
using Batsrus, PyPlot

file = "y.out"
bd = load(file)

DN = matplotlib.colors.DivergingNorm
set_cmap("RdBu_r")

contourf(bd, "uxS0", 50, plotrange=[-3,3,-3,3], plotinterval=0.05, norm=DN(0))
colorbar()
streamplot(bd, "uxS0;uzS0", density=2.0, color="g", plotrange=[-3,3,-3,3])
xlabel("x"); ylabel("y"); title("Ux [km/s]")

contourf(bd,"uxS0", 50, plotinterval=0.05, norm=DN(0))
colorbar()
axis("scaled")
xlabel("x"); ylabel("y"); title("uxS0")
```

### Tracing

The built-in `streamplot` function in Matplotlib is not satisfactory for accurately tracing. Instead we recommend [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) for tracing fieldlines and streamlines.

An example of tracing in a 2D cut and plot the field lines over contour:

```julia
using Batsrus, PyPlot

file = "test/y=0_var_1_t00000000_n00000000.out"
bd = load(file)

bx = bd.w[:,:,5]
bz = bd.w[:,:,7]
x  = bd.x[:,1,1]
z  = bd.x[1,:,2]

seeds = select_seeds(x, z; nSeed=100) # randomly select the seeding points

for i = 1:size(seeds)[2]
   xs = seeds[1,i]
   zs = seeds[2,i]
   # Tracing in both direction. Check the document for more options.
   x1, z1 = trace2d_eul(bx, bz, xs, zs, x, z, ds=0.1, maxstep=1000, gridType="ndgrid")
   plot(x1,z1,"--")
end
axis("equal")
```

Currently the `select_seeds` function uses pseudo random number generator that produces the same seeds every time.
