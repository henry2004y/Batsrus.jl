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

See the [examples](docs/src/man/examples.md).

## Guides

See [here](docs/src/man/guide.md) for some development thoughts.

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details.

## Acknowledgments

* All the nice guys who share their codes
