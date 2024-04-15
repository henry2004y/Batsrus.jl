# Data manipulation.

"""
	cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(bd::BATLData, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], dir::String="x", sequence::Int=1)

   x, w = bd.x, bd.w
   var_ = findfirst(x->x==lowercase(var), lowercase.(bd.head.wnames))
   isempty(var_) && error("$(var) not found in header variables!")

   if dir == "x"
      dim, d1, d2 = 1, 2, 3
   elseif dir == "y"
      dim, d1, d2 = 2, 1, 3
   else
      dim, d1, d2 = 3, 1, 2
   end

   cut1 = selectdim(view(x,:,:,:,d1), dim, sequence)
   cut2 = selectdim(view(x,:,:,:,d2), dim, sequence)
   W = selectdim(view(w,:,:,:,var_), dim, sequence)

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

   if isinf(limits[1]) limits[1] = hx[1] end
   if isinf(limits[3]) limits[3] = hy[1] end
   if isinf(limits[2]) limits[2] = hx[end] end
   if isinf(limits[4]) limits[4] = hy[end] end

   xind = searchsortedfirst(hx, limits[1]):searchsortedlast(hx, limits[2])
   yind = searchsortedfirst(hy, limits[3]):searchsortedlast(hy, limits[4])

   CartesianIndices((xind, yind))
end

@inline function findindexes(x, y, z, limits)
   hx = @view x[:,1,1]
   hy = @view y[1,:,1]
   hz = @view z[1,1,:]

   if isinf(limits[1]) limits[1] = hx[1] end
   if isinf(limits[3]) limits[3] = hy[1] end
   if isinf(limits[5]) limits[5] = hz[1] end
   if isinf(limits[2]) limits[2] = hx[end] end
   if isinf(limits[4]) limits[4] = hy[end] end
   if isinf(limits[6]) limits[6] = hz[end] end

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
    getvars(bd::BATLData, var::AbstractString) -> Array

Return variable data from string `var`. This is also supported via direct indexing,

# Examples
```
bd["rho"]
```

See also: [`getvars`](@ref).
"""
function getvar(bd::BATLData, var::AbstractString)
   if var in keys(variables_predefined)
      w = variables_predefined[var](bd)
   else
      var_ = findfirst(x->x==lowercase(var), lowercase.(bd.head.wnames))
      isnothing(var_) && error("$var not found in file header variables!")
      w = selectdim(bd.w, bd.head.ndim+1, var_)
   end

   w
end

@inline @Base.propagate_inbounds Base.getindex(bd::BATLData, var::AbstractString) =
   getvar(bd, var)

"""
    getvars(bd::BATLData, Names::Vector) -> Dict

Return variables' data as a dictionary from string vector.
See also: [`getvar`](@ref).
"""
function getvars(bd::BATLData{U}, Names::Vector{T}) where {U, T<:AbstractString}
   dict = Dict{T, Array{U}}()
   for name in Names
      dict[name] = getvar(bd, name)
   end

   dict
end

"Construct vectors from scalar components."
function _fill_vector_from_scalars(bd::BATLData, vstr1, vstr2, vstr3)
   v1 = getvar(bd, vstr1)
   v2 = getvar(bd, vstr2)
   v3 = getvar(bd, vstr3)
   v = Array{eltype(v1), ndims(v1)+1}(undef, 3, size(v1)...)

   Rpost = CartesianIndices(size(v1))
   for Ipost in Rpost
      v[1,Ipost] = v1[Ipost]
      v[2,Ipost] = v2[Ipost]
      v[3,Ipost] = v3[Ipost]
   end

   v
end

# Define derived parameters
const variables_predefined = Dict(
   "B2" => (bd -> @. $getvar(bd, "Bx")^2 + $getvar(bd, "By")^2 + $getvar(bd, "Bz")^2),
   "E2" => (bd -> @. $getvar(bd, "Ex")^2 + $getvar(bd, "Ey")^2 + $getvar(bd, "Ez")^2),
   "U2" => (bd -> @. $getvar(bd, "Ux")^2 + $getvar(bd, "Uy")^2 + $getvar(bd, "Uz")^2),
   "Ue2" => (bd -> @. $getvar(bd, "uxS0")^2 + $getvar(bd, "uyS0")^2 + $getvar(bd, "uzS0")^2),
   "Ui2" => (bd -> @. $getvar(bd, "uxS1")^2 + $getvar(bd, "uyS1")^2 + $getvar(bd, "uzS1")^2),
   "Bmag" => (bd -> @. sqrt($getvar(bd, "B2"))),
   "Emag" => (bd -> @. sqrt($getvar(bd, "E2"))),
   "Umag" => (bd -> @. sqrt($getvar(bd, "U2"))),
   "B" => (bd -> _fill_vector_from_scalars(bd, "Bx", "By", "Bz")),
   "E" => (bd -> _fill_vector_from_scalars(bd, "Ex", "Ey", "Ez")),
   "U" => (bd -> _fill_vector_from_scalars(bd, "Ux", "Uy", "Uz")),
)