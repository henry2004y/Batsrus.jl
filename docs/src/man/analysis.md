# Data Analysis

Once the data is loaded, you can extract specific variables and perform various analyses.

## Data Extraction

- Checking variable range

```julia
get_var_range(bd, "rho")
```

- Raw variables

Note that the variable names for queries must be in lowercase!

```julia
ρ = getvar(bd, "rho")
bd["rho"]
```

- Extracting data at a given location

```julia
loc = Float32[0.0, 0.0] # The type determines the output type
d = interp1d(bd, "rho", loc)
```

- Extracting data along a given line

```julia
point1 = Float32[-10.0, -1.0]
point2 = Float32[10.0, 1.0]
w = interp1d(bd, "rho", point1, point2)
```

```julia
bd["rho"][X=-10 .. 10, Y=Near(0.0)]
```

- Extracting data using DimensionalData

We can also use [Selectors](https://rafaqz.github.io/DimensionalData.jl/stable/selectors) from DimensionalData for extracting data. Note that the Selectors need to be imported from Batsrus.jl; alternatively you can simply `using DimensionalData`.

```julia
bd["rho"][X=At(0.0), Y=At(0.0)]
bd["rho"][X=-10 .. 10, Y=Near(0.5)]
bd["rho"][X=-10 .. 10, Y=-0.5 .. 0.5]
```

## Derived variables

We provide utility methods `get_magnitude`, `get_magnitude2`, `get_anisotropy`, `get_current_density`, `get_convection_E`, `get_hall_E`, and `get_pe_E` for calculating derived quantities. For convenience and performance, these can be accessed directly using symbols as indices:

```julia
# Magnitude access
Bmag = bd[:b]   # Equivalent to get_magnitude(bd, :B)
B2   = bd[:b2]  # Equivalent to get_magnitude2(bd, :B)
Emag = bd[:e]   # Equivalent to get_magnitude(bd, :E)
Umag = bd[:u]   # Equivalent to get_magnitude(bd, :U)
jmag = bd[:j]   # Equivalent to get_current_density magnitude

# Pressure anisotropy access
paniso0 = bd[:anisotropy0] # Equivalent to get_anisotropy(bd, 0)
paniso1 = bd[:anisotropy1] # Equivalent to get_anisotropy(bd, 1)

# Current density components
jx = bd[:jx]
jy = bd[:jy]
jz = bd[:jz]
```

Note that for PLANETARY units, the computed current density is automatically scaled to ``\\mu A/m^2``.

These symbol-based accesses are optimized to be near zero-allocation and are the recommended way to retrieve diagnostic variables in performance-critical loops.

### Electric fields

Electric field from the plasma simulations can be analyzed using the generalized Ohm's law.

- **Convection electric field**: $\mathbf{E}_{conv} = -\mathbf{u}_i \times \mathbf{B}$
- **Hall electric field**: $\mathbf{E}_{Hall} = (\mathbf{u}_i - \mathbf{u}_e) \times \mathbf{B}$
- **Electron pressure gradient electric field**: $\mathbf{E}_{p_e} = -\frac{1}{n_e e} \nabla \cdot \mathbf{P}_e$

```julia
# Convection E
Ecx, Ecy, Ecz = get_convection_E(bd)

# Hall E
Ehx, Ehy, Ehz = get_hall_E(bd)

# Electron pressure gradient E
# species=0 (default) is for electrons
Epx, Epy, Epz = get_pe_E(bd, species=0)
```

The `get_pe_E` function calculates the divergence of the electron pressure tensor. For structured grids, it uses central finite differences. The species mass is automatically retrieved from simulation parameters (e.g., `ms0`, `mass0`), but can be overridden with the `mass` keyword argument. Note that for `PLANETARY` or `km` units, the resulting electric field is scaled to $\mu V/m$ for consistency with the convection and Hall terms.

For structured grids, the current density $\mathbf{J} = \nabla \times \mathbf{B}$ can also be calculated using finite differences with `get_current_density`, which returns a tuple of `(Jx, Jy, Jz)`.

```julia
Jx, Jy, Jz = get_current_density(bd)
```

The underlying vector components can be retrieved using `get_vectors` or `get_vectors_indices`. `fill_vector_from_scalars` is also available for constructing a single array containing all three components:

```julia
Bvec = Batsrus.fill_vector_from_scalars(bd, :B)
```

Note that `fill_vector_from_scalars` involves additional array allocations compared to magnitude or anisotropy extraction.

Here is a full list of predefined derived quantities available via symbol indexing and `get_vectors`:

| Symbol / Derived Var | Meaning                          | Required variables |
|----------------------|----------------------------------|--------------------|
| :B                   | Magnetic field vector            | Bx, By, Bz         |
| :E                   | Electric field vector            | Ex, Ey, Ez         |
| :U                   | Velocity vector                  | Ux, Uy, Uz         |
| :U0                  | Species 0 velocity vector        | UxS0, UyS0, UzS0   |
| :U1                  | Species 1 velocity vector        | UxS1, UyS1, UzS1   |
| :J                   | Current density vector ($\mu A/m^2$)* | Bx, By, Bz         |
| :b, :b2              | Magnetic field magnitude (sq)    | Bx, By, Bz         |
| :e                   | Electric field magnitude         | Ex, Ey, Ez         |
| :u                   | Velocity magnitude               | Ux, Uy, Uz         |
| :j                   | Current density magnitude ($\mu A/m^2$)* | Bx, By, Bz         |
| :jx, :jy, :jz        | Current density components ($\mu A/m^2$)* | Bx, By, Bz         |
| :anisotropy0, 1, ... | Pressure anisotropy              | B and P tensor     |

\* Units are for PLANETARY data. For other systems, the units depend on the simulation normalization.

#### Finding indexes

To get the index of a certain quantity, e.g. electron number density

```julia
id = findindex(bd, "rhos0")
```

#### Get grid range

To extract the coordinates of the simulation grid:

```julia
ranges = get_range(bd)
# For 2D, this is equivalent to:
xrange, yrange = get_range(bd)
```

#### Get variable range

```julia
wmin, wmax = get_var_range(bd, var)
```
