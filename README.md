# VisAnaJulia
SWMF data reader and visualization using Julia.

This is inherited from the Matlab version of BATSRUS output reader and analyzer. It can be combined with the VTK format converter to generate files for Paraview and Tecplot.

The ultimate goal is to replace the IDL scripts for regular data visualizations, especially on Frontera. I am half way through. 

### Prerequisites

Julia 1.0+


## Usage
```
using Pkg; Pkg.activate(".")
using VisAna
```

IDL output processing:
- 1D binary 
```
filename = "1d_bin.out";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"p",plotmode="line")
plotdata(data[1],filehead[1],"p",plotmode="linegrid")
```

- 2D Cartesian (structured)
```
filename = "z=0_raw_1_t10.07620_n00000102.out";
filehead, data, filelist = readdata(filename,verbose=false);

plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar streamover")
plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar quiverover")
plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar streamover", density=2.0)
plotdata(data[1],filehead[1],"p",plotmode="grid")
plotdata(data[1],filehead[1],"p",plotmode="contbar",plotrange=[-50., 50., -1., 1.])
plotdata(data[1],filehead[1],"p",plotmode="contbar")
plotdata(data[1],filehead[1],"p",plotmode="contbarlog")
plotdata(data[1],filehead[1],"p",plotmode="surfbar")
```

- 2D unstructured
```
filename = "y=0_unstructured.outs";
filehead, data, filelist = readdata(filename, npict=2, verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar")
plotdata(data[1],filehead[1],"rho",plotmode="trimesh")
plotdata(data[1],filehead[1],"rho",plotmode="tricont")
```

- 2D structured spherical coordinates
```
filename = "y_structured.out"; #???
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar")
```

- 3D box
```
filename = "box.outs";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"bx",plotmode="contbar",cut="y")
plotdata(data[1],filehead[1],"bx",plotmode="contbar",cut="y",
   plotrange=[-1.4,-1.1,0.70,0.78])
```

- 3D structured spherical coordinates
```
filename = "3d_structured.out";
filehead, data, filelist = readdata(filename,verbose=false);
```

- log file
```
logfilename = "shocktube.log";
filehead, data = readlogdata(logfilename)
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

## Tricks

- This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own read! function in addition to the base ones.
- Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function `swaprows` is implemented to switch the orderings.

## Issues

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots.

- [ ] Switch to Makie for 3D plotting and animation
- [ ] PyBase support for manipulating data directly in Python

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details.

## Acknowledgments

* All the nice guys who share their codes
