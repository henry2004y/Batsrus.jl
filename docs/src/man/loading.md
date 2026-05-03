# Loading Data

Batsrus.jl supports various SWMF output formats, including IDL, log, and HDF5 files.

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

## HDF format output loader

```julia
using HDF5  # required to activate the HDF5 extension
filename = "3d__var_1_n00006288.h5"
file = BatsrusHDF5Uniform(filename)
```

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

## AMR Tree data

To load the block-adaptive tree structure for AMR datasets:

```julia
filetag = "3d__mhd_1_t00000000_n00000000"
batl = Batl(readhead(filetag * ".info"), readtree(filetag)...)
```

