using LinearAlgebra: tr

# Data manipulation.

"""
    cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(bd::BatsrusIDL, var::AbstractString;
      plotrange = [-Inf, Inf, -Inf, Inf], dir::String = "x", sequence::Int = 1)
   var_ = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname)
   isempty(var_) && error("$(var) not found in header variables!")

   if dir == "x"
      dim, d1, d2 = 1, 2, 3
   elseif dir == "y"
      dim, d1, d2 = 2, 1, 3
   else
      dim, d1, d2 = 3, 1, 2
   end

   cut1 = selectdim(view(bd.x, :, :, :, d1), dim, sequence)
   cut2 = selectdim(view(bd.x, :, :, :, d2), dim, sequence)
   W = selectdim(view(bd.w, :, :, :, var_), dim, sequence)

   if !all(isinf.(plotrange))
      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
   end

   cut1, cut2, W
end

@inline function checkvalidlimits(limits, dim::Int = 2)
   if dim == 2
      if length(limits) != 4
         throw(ArgumentError("Reduction range $limits should be [xmin xmax ymin ymax]!"))
      end

      if limits[1] > limits[2] || limits[3] > limits[4]
         throw(DomainError(limits, "Invalid reduction range!"))
      end
   elseif dim == 3
      if length(limits) != 6
         throw(
            ArgumentError(
            "Reduction range $limits should be [xmin xmax ymin ymax zmin max]!"),
         )
      end

      if limits[1] > limits[2] || limits[3] > limits[4] || limits[5] > limits[6]
         throw(DomainError(limits, "Invalid reduction range!"))
      end
   end
end

"""
    subsurface(x, y, data, limits)
    subsurface(x, y, u, v, limits)

Extract subset of 2D surface dataset in ndgrid format. See also: [`subvolume`](@ref).
"""
function subsurface(x, y, data, limits)
   checkvalidlimits(limits)

   if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
      DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup
      error("Selectors are only supported for rectilinear grids (gencoord=false)!")
   end

   xmin, xmax, ymin, ymax = limits
   xmin = isinf(xmin) ? first(dims(data, 1)) : xmin
   xmax = isinf(xmax) ? last(dims(data, 1)) : xmax
   ymin = isinf(ymin) ? first(dims(data, 2)) : ymin
   ymax = isinf(ymax) ? last(dims(data, 2)) : ymax

   selectors = (Between(xmin, xmax), Between(ymin, ymax))
   # This assumes x and y are the dimensions of data, which they should be for cutdata results
   subdata = data[selectors...]
   # subx and suby should also be sliced if they are DimArrays, or we can just return the dims from subdata
   # In cutdata, cut1 (x) and cut2 (y) are DimArrays.
   subx = x[selectors...]
   suby = y[selectors...]

   subx, suby, subdata
end

function subsurface(x, y, u, v, limits)
   checkvalidlimits(limits)

   if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
      DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup
      error("Selectors are only supported for rectilinear grids (gencoord=false)!")
   end

   xmin, xmax, ymin, ymax = limits
   xmin = isinf(xmin) ? first(dims(u, 1)) : xmin
   xmax = isinf(xmax) ? last(dims(u, 1)) : xmax
   ymin = isinf(ymin) ? first(dims(u, 2)) : ymin
   ymax = isinf(ymax) ? last(dims(u, 2)) : ymax

   selectors = (Between(xmin, xmax), Between(ymin, ymax))
   newu = u[selectors...]
   newv = v[selectors...]

   subx = x[selectors...]
   suby = y[selectors...]

   subx, suby, newu, newv
end

"""
    subvolume(x, y, z, data, limits)
    subvolume(x, y, z, u, v, w, limits)

Extract subset of 3D dataset in ndgrid format. See also: [`subsurface`](@ref).
"""
function subvolume(x, y, z, data, limits)
   checkvalidlimits(limits, 3)

   if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
      DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup ||
      DimensionalData.lookup(x, 3) isa DimensionalData.NoLookup
      error("Selectors are only supported for rectilinear grids (gencoord=false)!")
   end

   xmin, xmax, ymin, ymax, zmin, zmax = limits
   xmin = isinf(xmin) ? first(dims(data, 1)) : xmin
   xmax = isinf(xmax) ? last(dims(data, 1)) : xmax
   ymin = isinf(ymin) ? first(dims(data, 2)) : ymin
   ymax = isinf(ymax) ? last(dims(data, 2)) : ymax
   zmin = isinf(zmin) ? first(dims(data, 3)) : zmin
   zmax = isinf(zmax) ? last(dims(data, 3)) : zmax

   selectors = (Between(xmin, xmax), Between(ymin, ymax), Between(zmin, zmax))
   subdata = data[selectors...]
   subx = x[selectors...]
   suby = y[selectors...]
   subz = z[selectors...]

   subx, suby, subz, subdata
end

function subvolume(x, y, z, u, v, w, limits)
   checkvalidlimits(limits, 3)

   if DimensionalData.lookup(x, 1) isa DimensionalData.NoLookup ||
      DimensionalData.lookup(x, 2) isa DimensionalData.NoLookup ||
      DimensionalData.lookup(x, 3) isa DimensionalData.NoLookup
      error("Selectors are only supported for rectilinear grids (gencoord=false)!")
   end

   xmin, xmax, ymin, ymax, zmin, zmax = limits
   xmin = isinf(xmin) ? first(dims(u, 1)) : xmin
   xmax = isinf(xmax) ? last(dims(u, 1)) : xmax
   ymin = isinf(ymin) ? first(dims(u, 2)) : ymin
   ymax = isinf(ymax) ? last(dims(u, 2)) : ymax
   zmin = isinf(zmin) ? first(dims(u, 3)) : zmin
   zmax = isinf(zmax) ? last(dims(u, 3)) : zmax

   selectors = (Between(xmin, xmax), Between(ymin, ymax), Between(zmin, zmax))
   newu = u[selectors...]
   newv = v[selectors...]
   neww = w[selectors...]

   subx = x[selectors...]
   suby = y[selectors...]
   subz = z[selectors...]

   subx, suby, subz, newu, newv, neww
end

"""
     getvar(bd::BATS, var::AbstractString) -> Array

Return variable data from string `var`. This is also supported via direct indexing,

# Examples

```
bd["rho"]
```
"""
function getvar(
      bd::BatsrusIDL{ndim, TV}, var::AbstractString) where {ndim, TV}
   w = @view bd.w[var = At(lowercase(var))]
end

@inline Base.@propagate_inbounds Base.getindex(
bd::BatsrusIDL, var::AbstractString) = getvar(bd, var)

"""
     fill_vector_from_scalars(bd::BATS, var)

Construct vector of `var` from its scalar components. Alternatively, check
[`get_vectors`](@ref) for returning vector components as saparate arrays.
"""
function fill_vector_from_scalars(bd::BatsrusIDL, var)
   vt = get_vectors(bd, var)
   Rpost = CartesianIndices(size(bd.x)[1:(end - 1)])
   v = [vt[iv][i] for iv in 1:3, i in Rpost]
end

"""
     get_magnitude2(bd::BATS, var)

Calculate the magnitude square of vector `var`. See [`get_vectors`](@ref) for the options.
"""
function get_magnitude2(bd::BatsrusIDL, var = :B)
   vx, vy, vz = get_vectors(bd, var)
   v = similar(vx)

   @inbounds @simd for i in eachindex(v)
      v[i] = vx[i]^2 + vy[i]^2 + vz[i]^2
   end

   v
end

"""
     get_magnitude(bd::BATS, var)

Calculate the magnitude of vector `var`. See [`get_vectors`](@ref) for the options.
"""
function get_magnitude(bd::BatsrusIDL, var = :B)
   vx, vy, vz = get_vectors(bd, var)
   v = similar(vx)

   @inbounds @simd for i in eachindex(v)
      v[i] = √(vx[i]^2 + vy[i]^2 + vz[i]^2)
   end

   v
end

"""
     get_vectors(bd::BATS, var)

Return a tuple of vectors of `var`. `var` can be `:B`, `:E`, `:U`, or any `:U` followed by an index (e.g. `:U0` for species 0, `:U1` for species 1, etc.).
"""
function get_vectors(bd::BatsrusIDL, var)
   str = string(var)
   if str == "B"
      vx, vy, vz = bd["Bx"], bd["By"], bd["Bz"]
   elseif str == "E"
      vx, vy, vz = bd["Ex"], bd["Ey"], bd["Ez"]
   elseif str == "U"
      vx, vy, vz = bd["Ux"], bd["Uy"], bd["Uz"]
   else
      m = match(r"^U(\d+)$", str)
      if !isnothing(m)
         suffix = m[1]
         vx, vy, vz = bd["UxS" * suffix], bd["UyS" * suffix], bd["UzS" * suffix]
      else
         throw(ArgumentError("Vector variable $var not supported"))
      end
   end

   vx, vy, vz
end

"""
     get_anisotropy(bd::BATS, species=0)

Calculate the pressure anisotropy for `species`, indexing from 0. The default `method` is
based on the fact that the trace of the pressure tensor is a constant. The `rotation`
method is based on rotating the tensor.
"""
function get_anisotropy(bd::BatsrusIDL{2, TV}, species = 0;
      method::Symbol = :simple) where {TV}
   Bx, By, Bz = bd["Bx"], bd["By"], bd["Bz"]
   # Rotate the pressure tensor to align the 3rd direction with B
   pop = string(species)
   Pxx = bd["pXXS" * pop]
   Pyy = bd["pYYS" * pop]
   Pzz = bd["pZZS" * pop]
   Pxy = bd["pXYS" * pop]
   Pxz = bd["pXZS" * pop]
   Pyz = bd["pYZS" * pop]
   Paniso = similar(Pxx)

   @inbounds for j in axes(Pxx, 2), i in axes(Pxx, 1)
      b̂ = normalize(SA[Bx[i, j], By[i, j], Bz[i, j]])
      P = @SMatrix [Pxx[i, j] Pxy[i, j] Pxz[i, j];
                    Pxy[i, j] Pyy[i, j] Pyz[i, j];
                    Pxz[i, j] Pyz[i, j] Pzz[i, j]]

      if method == :simple
         p_parallel = b̂' * P * b̂
         p_perp = (tr(P) - p_parallel) / 2
         Paniso[i, j] = p_perp / p_parallel
      elseif method == :rotation
         Prot = rotateTensorToVectorZ(P, b̂)
         Paniso[i, j] = (Prot[1, 1] + Prot[2, 2]) / (2 * Prot[3, 3])
      else
         error("Unknown method for get_anisotropy: $method. Use :simple or :rotation.")
      end
   end

   Paniso
end

"""
Return the convection electric field from PIC outputs.
"""
function get_convection_E(bd::BatsrusIDL)
   Bx, By, Bz = get_vectors(bd, :B)
   # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
   uix, uiy, uiz = get_vectors(bd, :U1)

   Econvx = similar(Bx)
   Econvy = similar(By)
   Econvz = similar(Bz)
   # -Ui × B
   @simd for i in eachindex(Econvx)
      Econvx[i] = -uiy[i] * Bz[i] + uiz[i] * By[i]
      Econvy[i] = -uiz[i] * Bx[i] + uix[i] * Bz[i]
      Econvz[i] = -uix[i] * By[i] + uiy[i] * Bx[i]
   end

   Econvx, Econvy, Econvz
end

"""
Return the Hall electric field from PIC outputs.
"""
function get_hall_E(bd::BatsrusIDL)
   Bx, By, Bz = get_vectors(bd, :B)
   uex, uey, uez = get_vectors(bd, :U0)
   # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
   uix, uiy, uiz = get_vectors(bd, :U1)

   Ehallx = similar(Bx)
   Ehally = similar(By)
   Ehallz = similar(Bz)
   # (Ui - Ue) × B
   for i in eachindex(Ehallx)
      Ehallx[i] = (uiy[i] - uey[i]) * Bz[i] - (uiz[i] - uez[i]) * By[i]
      Ehally[i] = (uiz[i] - uez[i]) * Bx[i] - (uix[i] - uex[i]) * Bz[i]
      Ehallz[i] = (uix[i] - uex[i]) * By[i] - (uiy[i] - uey[i]) * Bx[i]
   end

   Ehallx, Ehally, Ehallz
end

"""
     get_timeseries(files::AbstractArray, loc; tstep = 1.0)

Extract plasma moments and EM field from PIC output `files` at `loc` with nearest neighbor.
Currently only works for 2D outputs. If a single point variable is needed, see [`interp1d`](@ref).
"""
function get_timeseries(files::AbstractArray, loc; tstep = 1.0)
   nfiles = length(files)
   bd = files[1] |> Batsrus.load
   xrange, yrange = get_range(bd)
   trange = range(bd.head.time, step = tstep, length = nfiles)
   @assert xrange[1]≤loc[1]≤xrange[end] "x location out of range!"
   @assert yrange[1]≤loc[2]≤yrange[end] "y location out of range!"
   x_ = searchsortedfirst(xrange, loc[1])
   y_ = searchsortedfirst(yrange, loc[2])
   v = zeros(Float32, 20, nfiles)

   @showprogress dt=1 desc="Extracting..." for it in eachindex(files)
      bd = files[it] |> Batsrus.load
      v[1, it] = bd["rhos0"][x_, y_]
      v[2, it] = bd["rhos1"][x_, y_]
      v[3, it] = bd["uxs0"][x_, y_]
      v[4, it] = bd["uys0"][x_, y_]
      v[5, it] = bd["uzs0"][x_, y_]
      v[6, it] = bd["uxs1"][x_, y_]
      v[7, it] = bd["uys1"][x_, y_]
      v[8, it] = bd["uzs1"][x_, y_]
      v[9, it] = bd["pxxs0"][x_, y_]
      v[10, it] = bd["pyys0"][x_, y_]
      v[11, it] = bd["pzzs0"][x_, y_]
      v[12, it] = bd["pxxs1"][x_, y_]
      v[13, it] = bd["pyys1"][x_, y_]
      v[14, it] = bd["pzzs1"][x_, y_]
      v[15, it] = bd["Bx"][x_, y_]
      v[16, it] = bd["By"][x_, y_]
      v[17, it] = bd["Bz"][x_, y_]
      v[18, it] = bd["Ex"][x_, y_]
      v[19, it] = bd["Ey"][x_, y_]
      v[20, it] = bd["Ez"][x_, y_]
   end

   trange, v
end

function get_timeseries(files::AbstractArray, loc::Vector; tstep = 1.0)
   get_timeseries(files, SVector{length(loc)}(loc); tstep)
end
