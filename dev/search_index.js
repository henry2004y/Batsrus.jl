var documenterSearchIndex = {"docs":
[{"location":"man/types/#Public-types-1","page":"Public types","title":"Public types","text":"","category":"section"},{"location":"man/types/#Public-types-in-module-Batsrus:-1","page":"Public types","title":"Public types in module Batsrus:","text":"","category":"section"},{"location":"man/types/#","page":"Public types","title":"Public types","text":"Modules = [Batsrus]\nPrivate = false\nOrder = [:type]","category":"page"},{"location":"man/types/#Batsrus.Data","page":"Public types","title":"Batsrus.Data","text":"Data\n\nPrimary data storage type, with fields head of header info, grid x, value w, and file info list.\n\n\n\n\n\n","category":"type"},{"location":"man/examples/#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"man/examples/#IDL-format-output-loader-1","page":"Examples","title":"IDL format output loader","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Read data","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"1d_bin.out\";\ndata = readdata(filename);\ndata = readdata(filename, verbose=true);\ndata = readdata(filename, npict=1);\ndata = readdata(filename, dir=\".\");","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"3D structured spherical coordinates","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_structured.out\";\ndata = readdata(filename, verbose=false);","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"log file","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"logfilename = \"shocktube.log\";\nhead, data = readlogdata(logfilename)","category":"page"},{"location":"man/examples/#Derived-variables-1","page":"Examples","title":"Derived variables","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"v = get_vars(data, [\"Bx\", \"By\", \"Bz\"])\nB = @. sqrt(v.Bx^2 + v.By^2 + v.Bz^2)","category":"page"},{"location":"man/examples/#Output-format-conversion-1","page":"Examples","title":"Output format conversion","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"We can convert 2D/3D BATSRUS outputs *.dat to VTK formats. The default converted filename is out.vtu.","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"ASCII Tecplot file (supports both tec and tcp):","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"x=0_mhd_1_n00000050.dat\"\nhead, data, connectivity = readtecdata(filename)\nconvertVTK(head, data, connectivity, outname)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_ascii.dat\"\nhead, data, connectivity = readtecdata(filename)\nconvertVTK(head, data, connectivity, outname)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Binary Tecplot file (set DOSAVETECBINARY=TRUE in BATSRUS PARAM.in):","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_bin.dat\"\nhead, data, connectivity = readtecdata(filename)\nconvertVTK(head, data, connectivity, outname)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"3D structured IDL file (gridType=1 returns rectilinear vtr file, gridType=2 returns structured vts file):","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_structured.out\"\nconvertBox2VTK(filename, gridType=1)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Multiple files:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"using Glob\nfilenamesIn = \"3d*.dat\"\ndir = \".\"\nfilenames = Vector{String}(undef,0)\nfilesfound = glob(filenamesIn, dir)\nfilenames = vcat(filenames, filesfound)\ntec = readtecdata.(filenames) # head, data, connectivity\nfor (i, outname) in enumerate(filenames)\n   convertVTK(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])\nend","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"If each individual file size is large, consider using:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"using Glob\nfilenamesIn = \"3d*.dat\"\ndir = \".\"\nfilenames = Vector{String}(undef,0)\nfilesfound = glob(filenamesIn, dir)\nfilenames = vcat(filenames, filesfound)\nfor (i, outname) in enumerate(filenames)\n   head, data, connectivity = readtecdata(outname)\n   convertVTK(head, data, connectivity, outname[1:end-4])\nend","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Multiple files in parallel:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"using Distributed\n@everywhere using Pkg\n@everywhere Pkg.activate(\"VisAnaJulia\");\n@everywhere using VisAna, Glob\n\nfilenamesIn = \"cut*.dat\"\ndir = \".\"\nfilenames = Vector{String}(undef,0)\nfilesfound = glob(filenamesIn, dir)\nfilenames = vcat(filenames, filesfound)\n\n@sync @distributed for outname in filenames\n   println(\"filename=$(outname)\")\n   head, data, connectivity = readtecdata(outname)\n   convertVTK(head, data, connectivity, outname[1:end-4])\nend","category":"page"},{"location":"#Batsrus.jl-Documentation-1","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"","category":"section"},{"location":"#Overview-1","page":"Batsrus.jl Documentation","title":"Overview","text":"","category":"section"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"note: Note\nThis package is still under development, so be careful for any future breaking changes!","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"BATSRUS data reader and converter in Julia.","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"This package inherits the ideas and code structures from its predecessor in IDL (developed by Gábor Tóth) and MATLAB (developed by Hongyang Zhou), and was originally part of VisAna. It uses the VTK XML format writer writeVTK to generate files for Paraview and Tecplot.","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"This package provides the following functionalities:","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"simulation data reader\ndata format conversion\nprogramming language interoperability","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"The ultimate goal is to build a convenient tool of reading and analyzing simulation outputs which is easy to install and easy to use.","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"tip: Ready to use?\nFeel free to contact the author for any help or collaboration!","category":"page"},{"location":"#Installation-1","page":"Batsrus.jl Documentation","title":"Installation","text":"","category":"section"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"Install VisAna from the julia REPL prompt with","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"using Pkg\nPkg.add(\"Batsrus\")","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"Pages = [\n    \"man/log.md\",\n    \"man/examples.md\",\n    \"man/functions.md\",\n    \"man/types.md\"\n]\nDepth = 1","category":"page"},{"location":"#Benchmark-1","page":"Batsrus.jl Documentation","title":"Benchmark","text":"","category":"section"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"Data loading speed of a 2.4GB 3D binary file, 317MB 3D binary file, and 65KB 2D binary file on Macbook Pro with quad core 2.2 GHz Intel i7 and 16 GB 1600 MHz DDR3:","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"2.4GB tmax [s] tmean [s]\nJulia 2.73 1.32\nPython 1.35 1.34\nIDL 6.18 6.08\nMATLAB 16.02 10.60","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"317MB tmean [ms]\nJulia 180.8\nPython 179.5\nIDL 453.5\nMATLAB 698.4","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"65KB tmean [μs]\nJulia 163.36\nPython 4390.95\nIDL 1970.29\nMATLAB 19273.25","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"The Julia, IDL, and MATLAB version all shares the same kernel design. The timings are obtained for Julia v1.3.1, Python 3.7.6 + Numpy 1.18.1, IDL 8.5, and MATLAB R2018b. For dynamic languages, the first time when function gets executed is usually also the slowest. Currently spacepy performs slightly better because of the well-optimized numpy library in C. For small data sizes, Julia is much faster than others.","category":"page"},{"location":"#Calling-From-Python-1","page":"Batsrus.jl Documentation","title":"Calling From Python","text":"","category":"section"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"In Python, you can easily take advantage of this package with the aid of PyJulia. After the installation, in the Python repl:","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"from julia import Batsrus\ndir = 'test'\nfilename = '1d__raw_2_t25.60000_n00000258.out'\ndata = Batsrus.readdata(filename, dir=dir)","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"There you have it! Enjoy!","category":"page"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"warning: Python dependency\nPyPlot package backend may be affected by the settings of PyJulia dependencies. If you want to set it back properly, you need to recompile the PyCall package in Julia.","category":"page"},{"location":"#Developers-1","page":"Batsrus.jl Documentation","title":"Developers","text":"","category":"section"},{"location":"#","page":"Batsrus.jl Documentation","title":"Batsrus.jl Documentation","text":"VisAna is developed by Hongyang Zhou.","category":"page"},{"location":"man/log/#Development-Log-1","page":"Development Log","title":"Development Log","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"All the workflows here is not restricted to one type of model output. After being familiar with new ideas and new models, one can easily make use of existing samples and create reader of their own. Because of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own read! function in addition to the base ones. Before v0.5.1, readdata function in Matlab for large data is 2 times faster than that in Julia. The reason is simply using read or unsafe_read in the lower level. The latter one is much faster. After the fix, Julia version performs 5 times faster than the Matlab version in reading binary data.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"I am enlightened by the way SpacePy handles files. Instead of a separate header and data array, it may be better to build a more contained struct. Also, the header could use NamedTuple instead of Dict.","category":"page"},{"location":"man/log/#Reading-Multiple-Files-1","page":"Development Log","title":"Reading Multiple Files","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"And actually, by far there is no single use case where I need to read multiple files together. If you want to do so, just call the function twice.","category":"page"},{"location":"man/log/#Interoperability-1","page":"Development Log","title":"Interoperability","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scatterred interpolation is done via Interpolate in Scipy. Hopefully these additional dependencies will be cut down.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"In the current version of PyCall and PyJulia, there is already direct support for accessing Julia struct objects (noted as jlwrap).","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"I have a new issue coming up with the interoperability with Python. I may need to split this package into pure IO and pure plotting to avoid the cross-dependency of Matplotlib. The idea is that PyPlot is only needed when I want to quickly scan through the data!","category":"page"},{"location":"man/log/#Units-1","page":"Development Log","title":"Units","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"I learned from Yuxi that the YT package can handle automatic unit conversion, just like what I saw from the Unitful.jl package.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"There is a unit package in Julia unitful for handling units. Take a look at that one if you really want to solve the unit problems.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"The ideas of how to make abstractions is more important than the unit conversion itself.","category":"page"},{"location":"man/log/#C-Dependency-1","page":"Development Log","title":"C Dependency","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"A real open-source project is a collaborated work not only from a bunch of people, but also a group of languages. In Julia, this can be achieved with the help of the Package manager.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"I want to have some C dependencies in my code instead of rewriting everything in Julia. This would serve as an attempt to quickly make things work.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Right now this seems to be a little bit difficult for me. I need to learn from experts. The tracing scheme in C is rewritten in Julia so I don't need to bother for now. Checkout BinaryBuilder for more information. A nice example is given in this C package.","category":"page"},{"location":"man/log/#VTK-1","page":"Development Log","title":"VTK","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"The VTK files does not have timestep information. To allow for further time series processing in Paraview, a script create_pvd.jl is provided for generating the pvd container. After dicussing with the author of WriteVTK.jl, this is finally achieved by adding a scalar that is not related to grid in VTK.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Currently, multi-block (VTM), rectilinear (VTR), and unstructured (VTU) conversion are supported. By default the file size will be reduced with compression level 6, but the actual compression ratio depends on the original data.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"For the plotting, streamline tracing and particle tracing, a common problem is the grid and related interpolation process. I am envisioning a more general approach to deal with block-based and unstructured grid to provide fundamental support for all of these.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Currently, BATSRUS output contains only cell center grid and cell center values. The multiblock VTK format conversion works, but it is incorrectly recognized as node center data. An ideal way is to output the original node coordinates and corresponding cell center values.","category":"page"},{"location":"man/log/#Ordering-of-connectivity-1","page":"Development Log","title":"Ordering of connectivity","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function swaprows is implemented to switch the orderings.","category":"page"},{"location":"man/log/#Variable-naming-1","page":"Development Log","title":"Variable naming","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Vector naming is messed up if you are using Tecplot VTK reader. For example, \"B [nT]\" –> \"B [nT]X\", \"B [nT]Y\", \"B [nT]_Z\". Not a big issue, but annoying.","category":"page"},{"location":"man/log/#AUXDATA-1","page":"Development Log","title":"AUXDATA","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"I have encountered a very bad problem of corrupting binary *.vtu files. It turned out that the issue is the starting position of data is wrong because of the way I skip the header AUXDATA part. Sometimes the binary numbers may contain newline character that confuses the reader. It is now fixed. Later on the reading of the header part is completely rewritten to provide better support for a variety of Tecplot Ascii headers. All the AXUDATA information is now stored into global VTK data.","category":"page"},{"location":"man/log/#Native-VTK-output-1","page":"Development Log","title":"Native VTK output","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"In the future versions of BATSRUS, we should be able to output VTK files directly with VTKFortran. I won't do it now.","category":"page"},{"location":"man/log/#Custom-VTK-reader-1","page":"Development Log","title":"Custom VTK reader","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"ParaView allows for custom Python reader. Examples can be found in Chapter 12.4 in the official manual, and an example of full Python plugins can be found at Kiware's gitlab page.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"The XML package not only provide writing into XML files, but also reading XML structures. Therefore, if you want you can create a VTK reader.","category":"page"},{"location":"man/log/#AMR-Grid-Structure-1","page":"Development Log","title":"AMR Grid Structure","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"vtkOverlappingAMR implements a somewhat strict Berger-Collela AMR scheme:","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"All grids are Cartesian.\nGrids at the same level do not overlap.\nThe refinement ratios, RL, between adjacent levels are integer (typically 2 or 4) and uniform within the same level.\nGrid cells are never partially refined; i.e., each cell is refined to four quads in 2D or eight hexahedra in 3D.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Or in other words,","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Refinement ratio across levels is constant.\nEach block at levels > 0 need to be covered 100% by one parent block of","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"previous level.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Some other restriction about what happens at the boundary.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"You can directly use vtkUniformGridAMR, which does not impose any restrictions. Most filters should work for this class - there just wouldn't be any specialized filters such as the dual-grid contour / clip ones for the vtkOverlappingAMR.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"The vtkAMRInformation documentation consists only of","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Refinement ratio between AMR levels\nGrid spacing for each level\nThe file block index for each block parent child information, if requested","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"(Image: sample_2DAMR) Sample 2D AMR Dataset with two levels and refinement ratio, RL=4. The root level (L0) consists of a single grid shown in black wireframe while the next level (L1) consists of two grids, depicted in green wireframe and red wireframe respectively. The two grids at L1 are projected from the root level to illustrate that the cells underneath are “hidden.”","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"In VTK, the collection of AMR grids is stored in a vtkHierarchicalBoxDataSet data-structure. Each grid, G(Li,k), is represented by a vtkUniformGrid data structure where the unique key pair (Li,k) denotes the corresponding level (Li) and the grid index within the level (k) with respect to the underlying hierarchical structure. An array historically known as IBLANK, stored as a cell attribute in vtkUniformGrid, denotes whether a cell is hidden or not. The blanking array is subsequently used by the mapper to hide lower resolution cells accordingly when visualizing the dataset.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"To enable the execution of data queries without loading the entire dataset in memory, metadata information is employed. The metadata stores a minimal set of geometric information for each grid in the AMR hierarchy. Specifically, the AMR metadata, B(Li,k), corresponding to the grid G(Li,k), is represented using a vtkAMRBox object and it consists of the following information:","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"N={Nx, Ny, Nz} — the cell dimensions of the grid (since the data is cell-centered)\nThe grid spacing at level L, hL={hx,hy,hz}\nThe grid level Li and grid index k\nThe global dataset origin, X=(X0, Y0, Z0), i.e., the minimum origin from all grids in level L0\nThe LoCorner and HiCorner, which describe the low and high corners of the rectangular region covered by the corresponding grid in a virtual integer lattice with the same spacing (h) that covers the entire domain.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"(Image: sample_2DAMR)","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Given the metadata information stored in the AMR box of each grid, the refinement ratio at each level can be easily computed using relationship (1) from Table 1. Further, the cartesian bounds the corresponding grid covers and the number of points and cells is also available (see relationships 2-4 in Table 1). Notably, geometric queries such as determining which cell contains a given point, or if a grid intersects a user-supplied slice plane, can be answered using just the metadata.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"There is a vtkAMRDualExtractionFilter, which constructs a dual-mesh (i.e., the mesh constructed by connecting the cell-centers) over the computational domain. If we can directly tell ParaView that the mesh we have is a dual-mesh, then the initial trial with multi-block data may work directly.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"AMRGaussianPulseSource","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"See Multi-Resolution Rendering with Overlapping AMR for the implementation of C++ reader in VTK.","category":"page"},{"location":"man/log/#Support-on-Derived-Variables-1","page":"Development Log","title":"Support on Derived Variables","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"Right now the derived quantity plots are not supported. In order to achieve this, I may need:","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"[x] A new function get_var(data, filehead, string) returning the derived variable\n[ ] A new plotting function that understands the derived data type","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct. This looks like the Python Calculator in ParaView.","category":"page"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"This test feature is dropped for now. ","category":"page"},{"location":"man/log/#Extra-Notes-1","page":"Development Log","title":"Extra Notes","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"I have already made a lot of mistakes by mixing the row-major and column-major codes. Explicitly list all the parts that require extra care!","category":"page"},{"location":"man/log/#Todo-List-1","page":"Development Log","title":"Todo List","text":"","category":"section"},{"location":"man/log/#","page":"Development Log","title":"Development Log","text":"[x] Full coverage of tests\n[x] PyJulia support for manipulating data directly in Python\n[x] Derived variable support (outdated)\n[x] General postprocessing script for concatenating and converting files.\n[x] Replace np.meshgrid with list comprehension\n[x] Drop the support for a long string containing several filenames; substitute by an array of strings.\n[ ] Find a substitution of triangulation in Julia\n[ ] Allow dot syntax to get dictionary contents (Base.convert?)\n[ ] Binary library support\n[ ] VTK reader","category":"page"},{"location":"man/functions/#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions/#Functions-exported-from-Batsrus:-1","page":"Functions","title":"Functions exported from Batsrus:","text":"","category":"section"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"Modules = [Batsrus]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"man/functions/#Batsrus.convertBox2VTK-Tuple{AbstractString}","page":"Functions","title":"Batsrus.convertBox2VTK","text":"convertBoxVTK(filename; dir=\".\", gridType=1, verbose=false)\n\nConvert 3D structured Tecplot data to VTK.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.convertVTK","page":"Functions","title":"Batsrus.convertVTK","text":"convertVTK(head, data, connectivity, filename=\"out\")\n\nConvert unstructured Tecplot data to VTK. Note that if using voxel type data in VTK, the connectivity sequence is different from Tecplot. Note that the 3D connectivity sequence in Tecplot is the same with the hexahedron type in VTK, but different with the voxel type. The 2D connectivity sequence is the same as the quad type, but different with the pixel type. For example, in 3D the index conversion is:\n\n# PLT to VTK voxel index_ = [1 2 4 3 5 6 8 7]\nfor i = 1:2\n   connectivity = swaprows!(connectivity, 4*i-1, 4*i)\nend\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#Batsrus.cutdata-Tuple{Data,AbstractString}","page":"Functions","title":"Batsrus.cutdata","text":"cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',\n\tcutPlaneIndex=1)\n\nGet 2D plane cut data of 3D box data.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.get_vars-Union{Tuple{T}, Tuple{Data,Array{T,1}}} where T<:AbstractString","page":"Functions","title":"Batsrus.get_vars","text":"Return data from input string vector.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.readdata-Tuple{AbstractString}","page":"Functions","title":"Batsrus.readdata","text":"readdata(filenameIn, (, dir=\".\", npict=1, verbose=false))\n\nRead data from BATSRUS output files. Stores the npict snapshot from an ascii or binary data file into the coordinates x and data w arrays. Filenames can be provided with wildcards.\n\nExamples\n\nfilename = \"1d_raw*\"\ndata = readdata(filename)\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.readlogdata-Tuple{AbstractString}","page":"Functions","title":"Batsrus.readlogdata","text":"Read information from log file.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.readtecdata-Tuple{AbstractString}","page":"Functions","title":"Batsrus.readtecdata","text":"readtecdata(filename, verbose=false)\n\nReturn header, data and connectivity from BATSRUS Tecplot outputs. Both 2D and 3D binary and ASCII formats are supported.\n\nExamples\n\nfilename = \"3d_ascii.dat\"\nhead, data, connectivity = readtecdata(filename)\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.showhead-Tuple{Batsrus.FileList,NamedTuple}","page":"Functions","title":"Batsrus.showhead","text":"showhead(file, filehead)\n\nDisplaying file header information.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.showhead-Tuple{Data}","page":"Functions","title":"Batsrus.showhead","text":"showhead(data)\n\nDisplaying file information for the Data type.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.subsurface-NTuple{4,Any}","page":"Functions","title":"Batsrus.subsurface","text":"subsurface(x, y, data, limits)\nsubsurface(x, y, u, v, limits)\n\nExtract subset of 2D surface dataset. This is a simplified version of subvolume.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#Batsrus.subvolume-NTuple{5,Any}","page":"Functions","title":"Batsrus.subvolume","text":"subvolume(x, y, z, data, limits)\nsubvolume(x, y, z, u, v, w, limits)\n\nExtract subset of 3D dataset in ndgrid format.\n\n\n\n\n\n","category":"method"}]
}
