# Data manipulation.

"""
	cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(bd::BATLData, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], dir::String="x", sequence::Int=1)
   var_ = findfirst(x->x==lowercase(var), lowercase.(bd.head.wnames))
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
    getvar(bd::BATLData, var::AbstractString) -> Array

Return variable data from string `var`. This is also supported via direct indexing,

# Examples
```
bd["rho"]
```
"""
function getvar(bd::BATLData{ndim, T}, var::AbstractString) where {ndim, T}
   if var in keys(variables_predefined)
      w = variables_predefined[var](bd)
   else
      var_ = findfirst(x->x==lowercase(var), lowercase.(bd.head.wnames))
      isnothing(var_) && error("$var not found in file header variables!")
      w = selectdim(bd.w, ndim+1, var_)
   end

   w
end

@inline @Base.propagate_inbounds Base.getindex(bd::BATLData, var::AbstractString) =
   getvar(bd, var)



"Construct vectors from scalar components."
function _fill_vector_from_scalars(bd::BATLData{ndim, T, U}, var) where {ndim, T, U}
   v1, v2, v3 = _get_vectors(bd, var)
   v = Array{T, ndims(v1)+1}(undef, 3, size(v1)...)

   Rpost = CartesianIndices(size(v1))
   for Ipost in Rpost
      v[1,Ipost] = v1[Ipost]
      v[2,Ipost] = v2[Ipost]
      v[3,Ipost] = v3[Ipost]
   end

   v
end

function _get_magnitude2(bd::BATLData{2, T, U}, var=:B) where {T, U}
   vx, vy, vz = _get_vectors(bd, var)
   v = similar(vx)::Array{T, 2}

   @simd for i in eachindex(v)
      v[i] = vx[i]^2 + vy[i]^2 + vz[i]^2
   end

   v
end

function _get_magnitude(bd::BATLData{2, T, U}, var=:B) where {T, U}
   vx, vy, vz = _get_vectors(bd, var)
   v = similar(vx)::Array{T, 2}

   @simd for i in eachindex(v)
      v[i] = √(vx[i]^2 + vy[i]^2 + vz[i]^2)
   end

   v
end

function _get_vectors(bd::BATLData, var)
   if var == :B
      vx = bd["Bx"]
      vy = bd["By"]
      vz = bd["Bz"]
   elseif var == :E
      vx = bd["Ex"]
      vy = bd["Ey"]
      vz = bd["Ez"]
   elseif var == :U
      vx = bd["Ux"]
      vy = bd["Uy"]
      vz = bd["Uz"]
   elseif var == :U0
      vx = bd["UxS0"]
      vy = bd["UyS0"]
      vz = bd["UzS0"]
   elseif var == :U1
      vx = bd["UxS1"]
      vy = bd["UyS1"]
      vz = bd["UzS1"] 
   end

   vx, vy, vz
end

function _get_anisotropy(bd::BATLData{2, T, U}, species=0) where {T, U}
   Bx, By, Bz = bd["Bx"], bd["By"], bd["Bz"]
   # Rotate the pressure tensor to align the 3rd direction with B
   pop = string(species)
   Pxx = bd["pXXS"*pop]
   Pyy = bd["pYYS"*pop]
   Pzz = bd["pZZS"*pop]
   Pxy = bd["pXYS"*pop]
   Pxz = bd["pXZS"*pop]
   Pyz = bd["pYZS"*pop]

   Paniso = similar(Pxx)::Array{T, 2}

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

# Define derived parameters
const variables_predefined = Dict{String, Function}(
   "B2" => (bd -> _get_magnitude2(bd, :B)),
   "E2" => (bd -> _get_magnitude2(bd, :E)),
   "U2" => (bd -> _get_magnitude2(bd, :U)),
   "Ue2" => (bd -> _get_magnitude2(bd, :U0)),
   "Ui2" => (bd -> _get_magnitude2(bd, :U1)),
   "Bmag" => (bd -> _get_magnitude(bd, :B)),
   "Emag" => (bd -> _get_magnitude(bd, :E)),
   "Umag" => (bd -> _get_magnitude(bd, :U)),
   "B" => (bd -> _fill_vector_from_scalars(bd, :B)),
   "E" => (bd -> _fill_vector_from_scalars(bd, :E)),
   "U" => (bd -> _fill_vector_from_scalars(bd, :U)),
   "Anisotropy0" => (bd -> _get_anisotropy(bd, 0)),
   "Anisotropy1" => (bd -> _get_anisotropy(bd, 1)),
)