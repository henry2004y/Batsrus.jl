# VisAnaJulia
SWMF data reader and visualization using Julia.

This is inherited from the Matlab version of BATSRUS output reader and analyzer. It can be combined with the VTK format converter to generate files for Paraview.

The ultimate goal is to replace the IDL scripts for regular data visualizations, especially on Frontera. I am half way through. 

### Prerequisites

Julia 1.0+

## Tricks

This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own read! function in addition to the base ones.

## Issues

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.m\
d) file for details

## Acknowledgments

* All the nice guys who share their codes
