# VisAna.jl Documentation

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

SWMF data reader and visualization using Julia.

This package is the inherited from its predecessor in IDL (developed by G.Toth) and Matlab (developed by H.Zhou).
It can be combined with the VTK format converter [writeVTK](https://github.com/jipolanco/WriteVTK.jl) to generate files for Paraview and Tecplot.
By default the file size will be reduced with compression level 6, but the actual compression ratio depends on the original data.

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scatterred interpolation is done via `Interpolate` in Scipy. Hopefully these additional dependencies will be cut down.

The ultimate goal is to build a convenient tool of reading and analyzing simulation outputs which is easy to install and easy to use.

!!! tip "Ready to use?"
    Feel free to contact the author for any help or collaboration!

## Installation
Install VisAna from the `julia` `REPL` prompt with
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/VisAnaJulia", rev="master"))
```

```@contents
Pages = [
    "man/guide.md",
    "man/examples.md",
    "man/functions.md",
    "man/types.md",
]
Depth = 1
```

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
