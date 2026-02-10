# APIs

## Batsrus

```@meta
CurrentModule = Batsrus
```

```@autodocs
Modules = [Batsrus]
```

## UnitfulBatsrus

```@meta
CurrentModule = Batsrus.UnitfulBatsrus
```

```@autodocs
Modules = [UnitfulBatsrus]
```

## HDF5 Extension

The HDF5 functionality is provided via a package extension. Load `HDF5` alongside `Batsrus` to enable it:

```julia
using Batsrus, HDF5
```

See [`BatsrusHDF5Uniform`](@ref) and [`extract_var`](@ref) for the API.
