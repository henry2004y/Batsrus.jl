# VisAna
SWMF data reader and visualization using Julia.

This is inherited from the Matlab version of BATSRUS output reader and analyzer. It can be combined with the VTK format converter to generate files for Paraview and Tecplot. By default the file size will be reduced with compression level 6, but the actual compression ratio depends on the original data.

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. For example, the 3D scatterred interpolation is done via `Interpolate` in Scipy.

The ultimate goal is to replace the IDL scripts for regular data visualizations, especially on Frontera. I am half way through.

### Prerequisites

Julia 1.0+

## Installation
```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/VisAnaJulia", rev="master"))
```

## Usage
```
#using Pkg; Pkg.activate(".") # for dev only
using VisAna
```

IDL format output processing:

- Read data
```
filename = "1d_bin.out";
filehead, data, filelist = readdata(filename, verbose=false);
filehead, data, filelist = readdata(filename);
filehead, data, filelist = readdata(filename, npict=1);
filehead, data, filelist = readdata(filename, dir=".");
```

A general `plotdata` function is provided for quick visualizations. In addition to that, some plotting functions can be directly called as shown below, which allows for more control by the user.

- 1D binary
```
plotdata(data[1], filehead[1], "p", plotmode="line")
plotdata(data[1], filehead[1], "p", plotmode="linegrid")
```

- 2D Cartesian (structured)
```
plotdata(data[1], filehead[1], "p bx;by", plotmode="contbar streamover")
plotdata(data[1], filehead[1], "p bx;by", plotmode="contbar quiverover")
plotdata(data[1], filehead[1], "p bx;by", plotmode="contbar streamover", density=2.0)
plotdata(data[1], filehead[1], "p", plotmode="grid")
plotdata(data[1], filehead[1], "p", plotmode="contbar", plotrange=[-50., 50., -1., 1.])
plotdata(data[1], filehead[1], "p", plotmode="contbar")
plotdata(data[1], filehead[1], "p", plotmode="contbarlog")
plotdata(data[1], filehead[1], "p", plotmode="surfbar")
```

- 2D unstructured
```
plotdata(data[1], filehead[1],"rho", plotmode="contbar")
plotdata(data[1], filehead[1],"rho", plotmode="trimesh")
plotdata(data[1], filehead[1],"rho", plotmode="tricont")
```

- 2D structured spherical coordinates
```
plotdata(data[1], filehead[1], "rho", plotmode="contbar")
```

- 3D box
```
plotdata(data[1], filehead[1], "bx", plotmode="contbar", cut="y")
plotdata(data[1], filehead[1], "bx", plotmode="contbar", cut="y", plotrange=[-1.4,-1.1,0.70,0.78])
```

- 3D structured spherical coordinates
```
filename = "3d_structured.out";
filehead, data, filelist = readdata(filename, verbose=false);
```

- log file
```
logfilename = "shocktube.log";
filehead, data = readlogdata(logfilename)
```

## Multiple dispatch for matplotlib functions
- line plot
```
c = plot(data[1], filehead[1], "p")
plt.setp(c, linestyle="--", linewidth=2);
```

- scatter plot
```
scatter(data[1], filehead[1], "p")
```

- contour
```
# 2D contour
contour(data[1], filehead[1], "p")
```

- filled contour
```
# 2D contourf
contourf(data[1], filehead[1], "p")
contourf(data[1], filehead[1], "p", levels, plotrange=[-10,10,-Inf,Inf], plotinterval=0.1)
```

- surface plot
```
# surface
plot_surface(data[1], filehead[1], "p")
```

- triangle surface plot
```
plot_trisurf(data[1], filehead[1], "p")
```

- triangle filled contour plot
```
tricontourf(data[1], filehead[1], "p")
```

## Derived variables
```
v = get_vars(data[1], filehead[1], ["Bx", "By", "Bz"])
B = @. sqrt(v.Bx^2 + v.By^2 + v.Bz^2)
```

## Output Format Conversion
ASCII tecplot file:
```
filename = "3d_ascii.dat"
head, data, connectivity  = readtecdata(filename,false)
convertVTK(head, data, connectivity, outname)
```

Binary tecplot file (`DOSAVETECBINARY=TRUE`):
```
filename = "3d_bin.dat"
head, data, connectivity  = readtecdata(filename,true)
convertVTK(head, data, connectivity, outname)
```

Multiple files:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
tec = readtecdata.(filenames, false) # head, data, connectivity
for (i, outname) in enumerate(filenames)
   convertVTK(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])
end
```

If each individual file size is large, consider using:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
for (i, outname) in enumerate(filenames)
   head, data, connectivity = readtecdata(filenames, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
```

Multiple files in parallel:
```
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("VisAnaJulia");
@everywhere using VisAna, Glob

filenamesIn = "cut*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
```

## Tricks

- This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own `read!` function in addition to the base ones.
- Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function `swaprows` is implemented to switch the orderings.
- Because of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.

## Issues

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots.

The current wrapper over Matplotlib makes it difficult to modify the plots afterwards, which especially causes problems when dealing with time series snapshots. The colorbar is so hard to fix.

In the roadmap of PyCall 2.0, there will direct support for accessing Julia objects. I hesitate to do it myself, so let's just wait for it to come.

The support for a long string containing several filenames as inputs has been dropped. It should be substituted by an array of strings.

Right now the derived quantity plots are not supported. In order to achieve this, I may need:
- [x] A new function `get_var(data, filehead, string)` returning the derived variable
- [ ] A new plotting function that understands the derived data type

The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct.

There is a user recipe in Plots. Check it out for the possibility of parameter control!

A direct wrapper over PyPlot function is possible, and would be more suitable for passing arguments. This may be a more plausible way to go than relying on recipes. Because of the reason described in [StackOverFlow](https://stackoverflow.com/questions/56153469/can-a-dictionary-be-used-as-a-complete-input-to-a-matplotlib-pyplot), plot cannot accept dictionary keywords input pairs.

When doing processing in batch mode on a cluster, there's usually no support for displaying backends. My current workaround:
```
using PyCall
matplotlib = pyimport("matplotlib")
matplotlib.use("Agg")
using PyPlot
```
This will give you a warning for the existed identifier, but it works.

Vector naming is messed up if you are using Tecplot VTK reader. For example, "B [nT]" --> "B [nT]_X", "B [nT]_Y", "B [nT]_Z". Not a big issue, but annoying.

- [ ] Fixed colorbar control through Matplotlib
- [ ] Test suite for checking validity
- [x] Cuts from 3D data visualization besides contour
- [ ] Switch to Makie for 3D plotting and animation
- [ ] PyBase support for manipulating data directly in Python
- [x] Derived variable support
- [ ] General postprocessing script for concatenating and converting files.
- [x] Direct wrapper over matplotlib functions to get seamless API
- [ ] Replace np.meshgrid with list comprehension

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details.

## Acknowledgments

* All the nice guys who share their codes
