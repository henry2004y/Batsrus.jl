# Batsrus.jl

[![Build Status](https://github.com/henry2004y/Batsrus.jl/workflows/CI/badge.svg)](https://github.com/henry2004y/Batsrus.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue.svg)][Batsrus-doc]
[![][codecov-img]][codecov-url]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4761843.svg)](https://doi.org/10.5281/zenodo.4761843)

Fast [BATSRUS](https://github.com/MSTEM-QUDA/BATSRUS)/[SWMF](https://github.com/MSTEM-QUDA/SWMF) data reading, converting, and visualizing using Julia, successor of [VisAnaMatlab](https://github.com/henry2004y/VisAnaMatlab).

This package provides the following functionalities:

* Simulation data reader
* 2D/3D region cut from the whole domain
* Interpolation from unstructured to structured data
* Data format conversion to VTK
* Data visualization via Makie, PyPlot, and Plots

For more details, please check the [document][Batsrus-doc].

## Prerequisites

* Julia 1.6+
* (Optional) Python and Matplotlib

## Installation

```julia
using Pkg
Pkg.add("Batsrus")
```

## Usage

See the [examples](https://henry2004y.github.io/Batsrus.jl/dev/man/examples/).

### Using from Python

In Python, you can easily take advantage of this package with the aid of [JuliaCall](https://juliapy.github.io/PythonCall.jl/dev/juliacall/) or [PyJulia](https://pyjulia.readthedocs.io/en/latest/).

With JuliaCall:

```python
from juliacall import Main as jl
jl.seval("using Batsrus")
file = 'example.out'
data = Batsrus.load(file, dir='test')
```

With PyJulia:

```python
from julia import Batsrus
file = 'example.out'
data = Batsrus.load(file, dir='test')
```

## Benchmark

See the [benchmark](https://henry2004y.github.io/Batsrus.jl/dev/#Benchmark-1) in the document.

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

[codecov-img]: https://codecov.io/gh/henry2004y/Batsrus.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/henry2004y/Batsrus.jl
[Batsrus-doc]: https://henry2004y.github.io/Batsrus.jl/dev
