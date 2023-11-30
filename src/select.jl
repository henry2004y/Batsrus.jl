# Data manipulation.

"""
	cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(data::BATLData, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], dir::String="x", sequence::Int=1)

   x, w = data.x, data.w
   var_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
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

"Return variable data from string `var`."
function getvar(data::BATLData, var::AbstractString)
   var_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
   isnothing(var_) && error("$var not found in file header variables!")

   w = selectdim(data.w, data.head.ndim+1, var_)
end

"""
    getvars(data::BATLData, Names::Vector) -> Dict

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


const variables_predefined = Dict(
   "B" => data -> sqrt.(getvar(data, "Bx").^2 .+ getvar(data, "By").^2 .+ getvar(data, "Bz").^2),
   "E" => data -> sqrt.(getvar(data, "Ex").^2 .+ getvar(data, "Ey").^2 .+ getvar(data, "Ez").^2),
   "U" => data -> sqrt.(getvar(data, "Ux").^2 .+ getvar(data, "Uy").^2 .+ getvar(data, "Uz").^2),
   #"beta" => data -> getvar(data, "P") ./ getvar(data, "B").^2 * 2μ,
)