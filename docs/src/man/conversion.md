# Output Format Conversion

We can convert 2D/3D BATSRUS outputs `*.dat` to VTK formats. It uses the VTK XML format writer [writeVTK](https://github.com/jipolanco/WriteVTK.jl) to generate files for Paraview and Tecplot. The default converted filename is `out.vtu`.

- ASCII Tecplot file (supports both `tec` and `tcp`) and binary Tecplot file (set `DOSAVETECBINARY=TRUE` in BATSRUS `PARAM.in`):

```julia
file = "x=0_mhd_1_n00000050.dat"
convertTECtoVTU(file)
```

- 3D structured IDL file (`gridType=:vti` returns image `vti` file, `gridType=:vtr` returns rectilinear `vtr` file, `gridType=:vts` returns structured `vts` file):

```julia
file = "3d_structured.out"
convertIDLtoVTK(file, gridType=:vti)
```

- 3D unstructured IDL file together with header and tree file:

```julia
filetag = "3d_var_1_n00002500"
convertIDLtoVTK(filetag)
```

!!! note
    The file suffix should not be provided for this to work correctly!

- Multiple files:

```julia
dir = "./"
filenames = filter(file -> startswith(file, "3d") && endswith(file, ".dat"), readdir(dir))
filenames = dir .* filenames

for filename in filenames
   convertTECtoVTU(filename, filename[1:end-4])
end
```

* Processing multiple files with threads in parallel:

```julia
dir = "./"
filenames = filter(file -> startswith(file, "3d") && endswith(file, ".dat"), readdir(dir))
filenames = dir .* filenames

Threads.@threads for filename in filenames
   println("filename=$filename")
   convertTECtoVTU(filename, filename[1:end-4])
end
```

More examples can be found in [examples](https://github.com/henry2004y/Batsrus.jl/tree/master/examples).
