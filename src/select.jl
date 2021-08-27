export getvars, getvar, cutdata, subvolume, subsurface


"""
	cutdata(data, var; plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

Get 2D plane cut in orientation `dir` for `var` out of 3D box `data` within `plotrange`.
The returned 2D data lies in the `sequence` plane from - to + in `dir`.
"""
function cutdata(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], dir="x", sequence=1)

   x, w = data.x, data.w
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

   X = @view x[:,:,:,1]
   Y = @view x[:,:,:,2]
   Z = @view x[:,:,:,3]

   W = w[:,:,:,VarIndex_]

   if dir == "x"
      cut1 = @view X[sequence,:,:]
      cut2 = @view Y[sequence,:,:]
      W    = @view W[sequence,:,:]
   elseif dir == "y"
      cut1 = @view X[:,sequence,:]
      cut2 = @view Z[:,sequence,:]
      W    = @view W[:,sequence,:]
   elseif dir == "z"
      cut1 = @view X[:,:,sequence]
      cut2 = @view Y[:,:,sequence]
      W    = @view W[:,:,sequence]
   end

   if !all(isinf.(plotrange))
      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
   end

   cut1, cut2, W
end

@inline function checkvalidlimits(limits, dim=2)
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

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])
   xind, yind
end

@inline function findindexes(x, y, z, limits)
   hx = @view x[:,1,1]
   hy = @view y[1,:,1]
   hz = @view z[1,1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[5]) limits[5] = minimum(hz) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end
   if isinf(limits[6]) limits[6] = maximum(hz) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])
   zind = findall(limits[5] .≤ hy .≤ limits[6])
   xind, yind, zind
end

"""
	subsurface(x, y, data, limits)
	subsurface(x, y, u, v, limits)

Extract subset of 2D surface dataset. See also: [`subvolume`](@ref).
"""
function subsurface(x, y, data, limits)

   checkvalidlimits(limits)

   xind, yind = findindexes(x, y, limits)

   newdata = subdata(data, xind, yind, size(data))

   newx = x[xind, yind]
   newy = y[xind, yind]

   newx, newy, newdata
end

function subsurface(x, y, u, v, limits)

   checkvalidlimits(limits)

   xind, yind = findindexes(x, y, limits)

   sz = size(u)
   newu = subdata(u, xind, yind, sz)
   newv = subdata(v, xind, yind, sz)

   newx = x[xind, yind]
   newy = y[xind, yind]

   newx, newy, newu, newv
end

"""
	subvolume(x, y, z, data, limits)
	subvolume(x, y, z, u, v, w, limits)

Extract subset of 3D dataset in ndgrid format. See also: [`subsurface`](@ref).
"""
function subvolume(x, y, z, data, limits)

   checkvalidlimits(limits, 3)

   xind, yind, zind = findindexes(x, y, z, limits)

   newdata = subdata(data, xind, yind, zind, size(data))

   newx = x[xind,yind,zind]
   newy = y[xind,yind,zind]
   newz = z[xind,yind,zind]

   newx, newy, newz, newdata
end

function subvolume(x, y, z, u, v, w, limits)

   checkvalidlimits(limits, 3)

   sz = size(u)

   xind, yind, zind = findindexes(x, y, z, limits)

   newu = subdata(u, xind, yind, zind, sz)
   newv = subdata(v, xind, yind, zind, sz)
   neww = subdata(w, xind, yind, zind, sz)

   newx = x[xind,yind,zind]
   newy = y[xind,yind,zind]
   newz = z[xind,yind,zind]

   newx, newy, newz, newu, newv, neww
end

"""
    subdata(data, xind, yind, sz)
    subdata(data, xind, yind, zind, sz)

Return the sliced data based on indexes `xind` and `yind` of size `sz`.
"""
function subdata(data, xind, yind, sz)
   newdata = data[xind,yind]
   newsz = size(newdata)

   if length(sz) > 2
      newdata = reshape(newdata, (newsz[1:2]..., sz[3:end]))
   end

   newdata
end

function subdata(data, xind, yind, zind, sz)

   newdata = data[xind,yind,zind]
   newsz = size(newdata)

   if length(sz) > 3
      newdata = reshape(newdata, (newsz[1:3]..., sz[4:end]))
   end

   newdata
end

"Return variable data from string `var`."
function getvar(data::Data, var)
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
   isnothing(VarIndex_) && error("$(var) not found in file header variables!")

   ndim = data.head.ndim
   if ndim == 1
      w = data.w[:,VarIndex_]
   elseif ndim == 2
      w = data.w[:,:,VarIndex_]
   elseif ndim == 3
      w = data.w[:,:,:,VarIndex_]
   end
   w
end

"""
    getvars(data::Data, Names::Vector) -> Dict

Return variables' data as a dictionary from string vector.
See also: [`getvar`](@ref).
"""
function getvars(data::Data, Names::Vector{T}) where T<:AbstractString

   dict = Dict()
   for name in Names
      dict[name] = getvar(data, name)
   end

   dict
end


const variables_predefined = Dict(
   "B" => data -> sqrt.(getvar(data, "Bx").^2 .+ getvar(data, "By").^2 .+ getvar(data, "Bz").^2),
   "E" => data -> sqrt.(getvar(data, "Ex").^2 .+ getvar(data, "Ey").^2 .+ getvar(data, "Ez").^2),
   "U" => data -> sqrt.(getvar(data, "Ux").^2 .+ getvar(data, "Uy").^2 .+ getvar(data, "Uz").^2),
   #"beta" => data -> getvar(data, "P") ./ getvar(data, "B").^2 * 2μ,
)