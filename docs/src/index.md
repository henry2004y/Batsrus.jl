# Batsrus.jl Documentation

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

BATSRUS data reader and converter in Julia.

This package inherits the ideas and code structures from its predecessor in IDL (developed by Gábor Tóth) and MATLAB (developed by Hongyang Zhou), and was originally part of [VisAna](https://github.com/henry2004y/VisAnaJulia).
It uses the VTK XML format writer [writeVTK](https://github.com/jipolanco/WriteVTK.jl) to generate files for Paraview and Tecplot.

This package provides the following functionalities:
  * simulation data reader
  * data format conversion
  * programming language interoperability

The ultimate goal is to build a convenient tool of reading and analyzing simulation outputs which is easy to install and easy to use.

!!! tip "Ready to use?"
    Feel free to contact the author for any help or collaboration!

## Installation
Install VisAna from the `julia REPL` prompt with
```julia
using Pkg
Pkg.add("Batsrus")
```

```@contents
Pages = [
    "man/log.md",
    "man/examples.md",
    "man/functions.md",
    "man/types.md"
]
Depth = 1
```

## Benchmark

Data loading speed of a 2.4GB 3D binary file, 317MB 3D binary file, and 65KB 2D binary file on Macbook Pro with quad core 2.2 GHz Intel i7 and 16 GB 1600 MHz DDR3:

| 2.4GB |   tmax [s] |  tmean [s] |
|:-------|:------:|:------:|
| Julia  | 2.73  |  1.32 |
| Python | 1.35  |  1.34 |
| IDL    | 6.18  |  6.08 |
| MATLAB | 16.02 | 10.60 |

| 317MB   | tmean [ms] |
|:-------|:---------:|
| Julia  | 180.8    |
| Python | 179.5   |
| IDL    | 453.5   |
| MATLAB | 698.4  |

| 65KB   | tmean [μs] |
|:-------|:---------:|
| Julia  | 163.36    |
| Python | 4390.95   |
| IDL    | 1970.29   |
| MATLAB | 19273.25  |

The Julia, IDL, and MATLAB version all shares the same kernel design. The timings are obtained for Julia v1.3.1, Python 3.7.6 + Numpy 1.18.1, IDL 8.5, and MATLAB R2018b.
For dynamic languages, the first time when function gets executed is usually also the slowest. Currently [spacepy](https://github.com/spacepy/spacepy) performs slightly better because of the well-optimized numpy library in C. For small data sizes, Julia is much faster than others.

## Calling From Python

In Python, you can easily take advantage of this package with the aid of [PyJulia](https://pyjulia.readthedocs.io/en/latest/).
After the installation, in the Python repl:
```python
from julia import Batsrus
dir = 'test'
filename = '1d__raw_2_t25.60000_n00000258.out'
data = Batsrus.readdata(filename, dir=dir)
```
There you have it! Enjoy!

!!! warning "Python dependency"
    `PyPlot` package backend may be affected by the settings of `PyJulia` dependencies. If you want to set it back properly, you need to recompile the `PyCall` package in Julia.

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
