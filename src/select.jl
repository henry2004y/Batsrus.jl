# Data manipulation.

"""
	cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(bd::BATS, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], dir::String="x", sequence::Int=1)
   var_ = findfirst(x->lowercase(x)==lowercase(var), bd.head.wname)
   isempty(var_) && error("$(var) not found in header variables!")

   if dir == "x"
      dim, d1, d2 = 1, 2, 3
   elseif dir == "y"
      dim, d1, d2 = 2, 1, 3
   else
      dim, d1, d2 = 3, 1, 2
   end

   cut1 = selectdim(view(bd.x,:,:,:,d1), dim, sequence)
   cut2 = selectdim(view(bd.x,:,:,:,d2), dim, sequence)
   W = selectdim(view(bd.w,:,:,:,var_), dim, sequence)

   if !all(isinf.(plotrange))
      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
   end

   cut1, cut2, W
end

@inline function checkvalidlimits(limits, dim::Int=2)
   if dim == 2
      if length(limits) != 4
         throw(ArgumentError("Reduction range $limits should be [xmin xmax ymin ymax]!"))
      end

      if limits[1] > limits[2] || limits[3] > limits[4]
         throw(DomainError(limits, "Invalid reduction range!"))
      end
   elseif dim == 3
      if length(limits) != 6
         throw(ArgumentError(
            "Reduction range $limits should be [xmin xmax ymin ymax zmin max]!"))
      end

      if limits[1] > limits[2] || limits[3] > limits[4] || limits[5] > limits[6]
         throw(DomainError(limits, "Invalid reduction range!"))
      end
   end
end

@inline function findindexes(x, y, limits)
   hx = @view x[:,1]
   hy = @view y[1,:]

   limits[1] = ifelse(isinf(limits[1]), hx[1], limits[1])
   limits[2] = ifelse(isinf(limits[2]), hx[end], limits[2])
   limits[3] = ifelse(isinf(limits[3]), hy[1], limits[3])
   limits[4] = ifelse(isinf(limits[4]), hy[end], limits[4])

   xind = searchsortedfirst(hx, limits[1]):searchsortedlast(hx, limits[2])
   yind = searchsortedfirst(hy, limits[3]):searchsortedlast(hy, limits[4])

   CartesianIndices((xind, yind))
end

@inline function findindexes(x, y, z, limits)
   hx = @view x[:,1,1]
   hy = @view y[1,:,1]
   hz = @view z[1,1,:]

   limits[1] = ifelse(isinf(limits[1]), hx[1], limits[1])
   limits[2] = ifelse(isinf(limits[2]), hx[end], limits[2])
   limits[3] = ifelse(isinf(limits[3]), hy[1], limits[3])
   limits[4] = ifelse(isinf(limits[4]), hy[end], limits[4])
   limits[5] = ifelse(isinf(limits[5]), hz[1], limits[5])
   limits[6] = ifelse(isinf(limits[6]), hz[end], limits[6])

   xind = searchsortedfirst(hx, limits[1]):searchsortedlast(hx, limits[2])
   yind = searchsortedfirst(hy, limits[3]):searchsortedlast(hy, limits[4])
   zind = searchsortedfirst(hz, limits[5]):searchsortedlast(hz, limits[6])

   CartesianIndices((xind, yind, zind))
end

"""
	subsurface(x, y, data, limits)
	subsurface(x, y, u, v, limits)

Extract subset of 2D surface dataset in ndgrid format. See also: [`subvolume`](@ref).
"""
function subsurface(x, y, data, limits)
   checkvalidlimits(limits)

   ids = findindexes(x, y, limits)

   @views begin
      subdata = data[ids]

      subx = x[ids]
      suby = y[ids]
   end

   subx, suby, subdata
end

function subsurface(x, y, u, v, limits)
   checkvalidlimits(limits)

   ids = findindexes(x, y, limits)

   @views begin
      newu = u[ids]
      newv = v[ids]

      subx = x[ids]
      suby = y[ids]
   end

   subx, suby, newu, newv
end

"""
	subvolume(x, y, z, data, limits)
	subvolume(x, y, z, u, v, w, limits)

Extract subset of 3D dataset in ndgrid format. See also: [`subsurface`](@ref).
"""
function subvolume(x, y, z, data, limits)
   checkvalidlimits(limits, 3)

   ids = findindexes(x, y, z, limits)

   @views begin
      subdata = data[ids]

      subx = x[ids]
      suby = y[ids]
      subz = z[ids]
   end

   subx, suby, subz, subdata
end

function subvolume(x, y, z, u, v, w, limits)
   checkvalidlimits(limits, 3)

   ids = findindexes(x, y, z, limits)

   @views begin
      newu = u[ids]
      newv = v[ids]
      neww = w[ids]

      subx = x[ids]
      suby = y[ids]
      subz = z[ids]
   end

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
function getvar(bd::BATS{ndim, TV, TX, TW}, var::AbstractString) where {ndim, TV, TX, TW}
   if var == "B2"
      w = get_magnitude2(bd, :B)
   elseif var == "E2"
      w = get_magnitude2(bd, :E)
   elseif var == "U2"
      w = get_magnitude2(bd, :U)
   elseif var == "Ue2"
      w = get_magnitude2(bd, :U0)
   elseif var == "Ui2"
      w = get_magnitude2(bd, :U1)
   elseif var == "Bmag"
      w = get_magnitude(bd, :B)
   elseif var == "Emag"
      w = get_magnitude(bd, :E)
   elseif var == "Umag"
      w = get_magnitude(bd, :U)
   elseif var == "Uemag"
      w = get_magnitude(bd, :U0)
   elseif var == "Uimag"
      w = get_magnitude(bd, :U1)
   elseif var == "B"
      w = fill_vector_from_scalars(bd, :B)
   elseif var == "E"
      w = fill_vector_from_scalars(bd, :E)
   elseif var == "U"
      w = fill_vector_from_scalars(bd, :U)
   elseif var == "Anisotropy0"
      w = get_anisotropy(bd, 0)
   elseif var == "Anisotropy1"
      w = get_anisotropy(bd, 1)
   else
      w = @view bd.w[var=At(lowercase(var))]
   end

   w
end

@inline @Base.propagate_inbounds Base.getindex(bd::BATS, var::AbstractString) =
   getvar(bd, var)


"""
    fill_vector_from_scalars(bd::BATS, var)

Construct vector of `var` from its scalar components. Alternatively, check
[`get_vectors`](@ref) for returning vector components as saparate arrays.
"""
function fill_vector_from_scalars(bd::BATS, var)
   vt = get_vectors(bd, var)
   Rpost = CartesianIndices(size(bd.x)[1:bd.head.ndim])
   v = [vt[iv][i] for iv in 1:3, i in Rpost]
end

"""
    get_magnitude2(bd::BATS, var)

Calculate the magnitude square of vector `var`. See [`get_vectors`](@ref) for the options.
"""
function get_magnitude2(bd::BATS, var=:B)
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
function get_magnitude(bd::BATS, var=:B)
   vx, vy, vz = get_vectors(bd, var)
   v = similar(vx)

   @inbounds @simd for i in eachindex(v)
      v[i] = √(vx[i]^2 + vy[i]^2 + vz[i]^2)
   end

   v
end

"""
    get_vectors(bd::BATS, var)

Return a tuple of vectors of `var`. `var` can be `:B`, `:E`, `:U`, `:U0`, or `:U1`.
"""
function get_vectors(bd::BATS, var)
   if var == :B
      vx, vy, vz = bd["Bx"], bd["By"], bd["Bz"]
   elseif var == :E
      vx, vy, vz = bd["Ex"], bd["Ey"], bd["Ez"]
   elseif var == :U
      vx, vy, vz = bd["Ux"], bd["Uy"], bd["Uz"]
   elseif var == :U0
      vx, vy, vz = bd["UxS0"], bd["UyS0"], bd["UzS0"]
   elseif var == :U1
      vx, vy, vz = bd["UxS1"], bd["UyS1"], bd["UzS1"]
   end

   vx, vy, vz
end

"""
    get_anisotropy(bd::BATS, species=0)

Calculate the pressure anisotropy for `species`, indexing from 0.
"""
function get_anisotropy(bd::BATS{2, TV, TX, TW}, species=0) where {TV, TX, TW}
   Bx, By, Bz = bd["Bx"], bd["By"], bd["Bz"]
   # Rotate the pressure tensor to align the 3rd direction with B
   pop = string(species)
   Pxx = bd["pXXS"*pop]
   Pyy = bd["pYYS"*pop]
   Pzz = bd["pZZS"*pop]
   Pxy = bd["pXYS"*pop]
   Pxz = bd["pXZS"*pop]
   Pyz = bd["pYZS"*pop]
   #TODO: Generalize to n-dimension with CartesianIndices
   Paniso = similar(Pxx)

   @inbounds for j in axes(Pxx, 2), i in axes(Pxx, 1)  
      b̂ = normalize(SA[Bx[i,j], By[i,j], Bz[i,j]])
      P =  @SMatrix [
         Pxx[i,j] Pxy[i,j] Pxz[i,j];
         Pxy[i,j] Pyy[i,j] Pyz[i,j];
         Pxz[i,j] Pyz[i,j] Pzz[i,j]]

      Prot = rotateTensorToVectorZ(P, b̂)
      Paniso[i,j] = (Prot[1,1] + Prot[2,2]) / (2*Prot[3,3])
   end

   Paniso
end

"Return the convection electric field from PIC outputs."
function get_convection_E(bd::BATS)
   Bx, By, Bz = get_vectors(bd, :B)
   # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
   uix, uiy, uiz = get_vectors(bd, :U1)

   Econvx = similar(Bx)
   Econvy = similar(By)
   Econvz = similar(Bz)
   # -Ui × B
   @simd for i in eachindex(Econvx)
      Econvx[i] = -uiy[i]*Bz[i] + uiz[i]*By[i]
      Econvy[i] = -uiz[i]*Bx[i] + uix[i]*Bz[i]
      Econvz[i] = -uix[i]*By[i] + uiy[i]*Bx[i]
   end

   Econvx, Econvy, Econvz
end

"Return the Hall electric field from PIC outputs."
function get_hall_E(bd::BATS)
   Bx, By, Bz = get_vectors(bd, :B)
   uex, uey, uez = get_vectors(bd, :U0)
   # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
   uix, uiy, uiz = get_vectors(bd, :U1)

   Ehallx = similar(Bx)
   Ehally = similar(By)
   Ehallz = similar(Bz)
   # (Ui - Ue) × B
   for i in eachindex(Ehallx)
      Ehallx[i] = (uiy[i] - uey[i])*Bz[i] - (uiz[i] - uez[i])*By[i]
      Ehally[i] = (uiz[i] - uez[i])*Bx[i] - (uix[i] - uex[i])*Bz[i]
      Ehallz[i] = (uix[i] - uex[i])*By[i] - (uiy[i] - uey[i])*Bx[i]
   end

   Ehallx, Ehally, Ehallz
end
