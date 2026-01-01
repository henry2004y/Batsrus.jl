# Utility functions for plotting and analyzing.

"""
     interp2d(bd::BATS, var::AbstractString, plotrange=[-Inf, Inf, -Inf, Inf],
    	 plotinterval=Inf; kwargs...)

Return 2D interpolated slices of data `var` from `bd`. If `plotrange` is not set, output
data resolution is the same as the original.

# Keyword Arguments

  - `innermask=false`: Whether to mask the inner boundary with NaN.
  - `rbody=1.0`: Radius of the inner mask. Used when the rbody parameter is not found in the header.
  - `useMatplotlib=true`: Whether to Matplotlib (faster) or NaturalNeighbours for scattered
    interpolation. If true, a linear interpolation is performed on a constructed triangle mesh.
"""
function interp2d(bd::BATS{2, TV, TX, TW}, var::AbstractString,
      plotrangeIn::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32;
      innermask::Bool = false, rbody::Real = 1.0, useMatplotlib::Bool = true
) where {TV, TX, TW}
   x, w = bd.x, bd.w
   varIndex_ = findindex(bd, var)
   plotrange = TV.(plotrangeIn)

   if bd.head.gencoord # Generalized coordinates
      X, Y = eachslice(x, dims = 3)
      X, Y = vec(X), vec(Y)
      W = @views w[:, :, varIndex_] |> vec

      adjust_plotrange!(plotrange, extrema(X), extrema(Y))
      # Set a heuristic value if not set
      if isinf(plotinterval)
         plotinterval = (plotrange[2] - plotrange[1]) / size(X, 1)
      end
      xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
      yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)

      if useMatplotlib
         try
            Wi = _triangulate_matplotlib(X, Y, W, xi, yi)
         catch e
            if e isa MethodError
               error("Matplotlib interpolation requires PyPlot to be loaded.")
            else
               rethrow(e)
            end
         end
      else
         xi, yi, Wi = interpolate2d_generalized_coords(X, Y, W, plotrange, plotinterval)
      end
   else # Cartesian coordinates
      xrange, yrange = get_range(bd)
      if all(isinf.(plotrange))
         xi, yi = xrange, yrange
         Wi = w[:, :, varIndex_].data' # Matplotlib does not accept view!
      else
         adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

         if isinf(plotinterval)
            xi = range(plotrange[1], stop = plotrange[2], step = xrange[2] - xrange[1])
            yi = range(plotrange[3], stop = plotrange[4], step = yrange[2] - yrange[1])
         else
            xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
            yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
         end
         itp = @views scale(interpolate(w[:, :, varIndex_], BSpline(Linear())),
            (xrange, yrange))
         Wi = [itp(i, j) for j in yi, i in xi]
      end
   end

   # Mask a circle at the inner boundary
   if innermask
      varIndex_ = findlast(x -> x == "rbody", bd.head.param)
      if isnothing(varIndex_)
         @info "rbody not found in file header parameters; use keyword rbody"
         @inbounds @simd for i in CartesianIndices(Wi)
            if xi[i[2]]^2 + yi[i[1]]^2 < rbody^2
               Wi[i] = NaN
            end
         end
      else
         ndim = 2
         ParamIndex_ = varIndex_ - ndim - bd.head.nw
         @inbounds @simd for i in CartesianIndices(Wi)
            if xi[i[1]]^2 + yi[i[2]]^2 < bd.head.eqpar[ParamIndex_]^2
               Wi[i] = NaN
            end
         end
      end
   end

   xi, yi, Wi
end

"""
Return the axis range for 2D outputs. See [`interp2d`](@ref).
"""
function meshgrid(bd::BATS,
      plotrange::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32)
   x = bd.x

   if bd.head.gencoord # Generalized coordinates
      X, Y = eachslice(x, dims = 3)
      X, Y = vec(X), vec(Y)

      adjust_plotrange!(plotrange, extrema(X), extrema(Y))
      # Set a heuristic value if not set
      if isinf(plotinterval)
         plotinterval = (plotrange[2] - plotrange[1]) / size(X, 1)
      end
      xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
      yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
   else # Cartesian coordinates
      xrange, yrange = get_range(bd)
      if all(isinf.(plotrange))
         xi, yi = xrange, yrange
      else
         adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

         if isinf(plotinterval)
            xi = range(plotrange[1], stop = plotrange[2], step = xrange[2] - xrange[1])
            yi = range(plotrange[3], stop = plotrange[4], step = yrange[2] - yrange[1])
         else
            xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
            yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
         end
      end
   end

   xi, yi
end

"""
Find variable index in the BATSRUS data.
"""
function findindex(bd::BATS, var::AbstractString)
   varIndex_ = findfirst(x -> lowercase(x) == lowercase(var), bd.head.wname)
   isnothing(varIndex_) && error("$(var) not found in file header variables!")

   varIndex_
end

"""
Generating consistent 2D arrays for passing to plotting functions.
"""
function meshgrid(x, y)
   X = [x for _ in y, x in x]
   Y = [y for y in y, _ in x]

   X, Y
end

@inline hasunit(bd::BATS) = startswith(bd.head.headline, "normalized") ? false : true

"""
Adjust 2D plot ranges.
"""
function adjust_plotrange!(plotrange, xlimit, ylimit)
   plotrange[1] = ifelse(isinf(plotrange[1]), xlimit[1], plotrange[1])
   plotrange[2] = ifelse(isinf(plotrange[2]), xlimit[2], plotrange[2])
   plotrange[3] = ifelse(isinf(plotrange[3]), ylimit[1], plotrange[3])
   plotrange[4] = ifelse(isinf(plotrange[4]), ylimit[2], plotrange[4])

   return
end

"""
Perform Triangle interpolation of 2D data `W` on grid `X`, `Y`.
"""
function interpolate2d_generalized_coords(X::T, Y::T, W::T,
      plotrange::Vector{<:AbstractFloat}, plotinterval::Real) where {T <: AbstractVector}
   xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
   yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
   itp = NN.interpolate(X, Y, W)
   #TODO Interpolate as a whole using predicates with multithreading
   Wi = [itp(x, y; method = NN.Triangle()) for y in yi, x in xi]::Matrix{eltype(W)}

   xi, yi, Wi
end

"""
     interp1d(bd::BATS, var::AbstractString, loc::AbstractVector{<:AbstractFloat})

Interpolate `var` at spatial point `loc` in `bd`.
"""
function interp1d(
      bd::BATS{2, TV, TX, TW},
      var::AbstractString,
      loc::AbstractVector{<:AbstractFloat}
) where {TV, TX, TW}
   @assert !bd.head.gencoord "Only accept structured grids!"

   v = getview(bd, var)
   xrange, yrange = get_range(bd)
   itp = scale(interpolate(v, BSpline(Linear())), (xrange, yrange))

   Wi = itp(loc...)
end

"""
     interp1d(bd::BATS, var::AbstractString, point1::Vector, point2::Vector)

Interpolate `var` along a line from `point1` to `point2` in `bd`.
"""
function interp1d(
      bd::BATS{2, TV, TX, TW},
      var::AbstractString,
      point1::Vector,
      point2::Vector
) where {TV, TX, TW}
   @assert !bd.head.gencoord "Only accept structured grids!"

   v = getview(bd, var)
   xrange, yrange = get_range(bd)
   itp = scale(interpolate(v, BSpline(Linear())), (xrange, yrange))
   lx = point2[1] - point1[1]
   ly = point2[2] - point1[2]
   nx = lx ÷ xrange.step |> Int
   ny = ly ÷ yrange.step |> Int
   ns = floor(Int, √(nx^2 + ny^2))
   dx = lx / ns
   dy = ly / ns
   points = [(point1[1] + i * dx, point1[2] + i * dy) for i in 0:ns]

   Wi = [itp(loc...) for loc in points]
end

"""
     slice1d(bd, var, icut::Int=1, dir::Int=2)

Return view of variable `var` in `bd` along 1D slice. `icut` is the index along axis `dir`.
`dir == 1` means align with the 2nd (y) axis, `dir == 2` means align with the 1st (x) axis.
"""
slice1d(bd, var, icut::Int = 1, dir::Int = 2) = selectdim(bd[var], dir, icut)

"""
Return view of variable `var` in `bd`.
"""
function getview(bd::BATS{1, TV, TX, TW}, var) where {TV, TX, TW}
   varIndex_ = findindex(bd, var)

   v = @view bd.w[:, varIndex_]
end

function getview(bd::BATS{2, TV, TX, TW}, var) where {TV, TX, TW}
   varIndex_ = findindex(bd, var)

   v = @view bd.w[:, :, varIndex_]
end

"""
Return value range of `var` in `bd`.
"""
get_var_range(bd::BATS, var) = getview(bd, var) |> extrema

"""
Return mesh range of `bd`.
"""
function get_range(bd::BATS{2, TV, TX, TW}) where {TV, TX, TW}
   x = bd.x
   xrange = range(x[1, 1, 1], x[end, 1, 1], length = size(x, 1))
   yrange = range(x[1, 1, 2], x[1, end, 2], length = size(x, 2))

   xrange, yrange
end

function get_range(bd::BATS{3, TV, TX, TW}) where {TV, TX, TW}
   x = bd.x
   xrange = range(x[1, 1, 1, 1], x[end, 1, 1, 1], length = size(x, 1))
   yrange = range(x[1, 1, 1, 2], x[1, end, 1, 2], length = size(x, 2))
   zrange = range(x[1, 1, 1, 3], x[1, 1, end, 3], length = size(x, 3))

   xrange, yrange, zrange
end

"""
Squeeze singleton dimensions for an array `A`.
"""
function squeeze(A::AbstractArray)
   singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)

   dropdims(A, dims = singleton_dims)
end

"""
     rotateTensorToVectorZ(tensor::AbstractMatrix, v::AbstractVector) -> SMatrix{3,3}

Rotate `tensor` with a rotation matrix that aligns the 3rd direction with `vector`, which is equivalent to change the basis from (i,j,k) to (i′,j′,k′) where k′ ∥ vector.
Reference: [Tensor rotation](https://math.stackexchange.com/questions/2303869/tensor-rotation)
"""
function rotateTensorToVectorZ(tensor::AbstractMatrix{T}, v) where T
   k = SVector{3, T}(0.0, 0.0, 1.0)
   axis = v × k |> normalize
   if axis[1] == axis[2] == 0
      return tensor
   else
      angle = acos(v ⋅ k / sqrt(v[1]^2 + v[2]^2 + v[3]^2))
      R = getRotationMatrix(axis, angle)
      return R * tensor * R'
   end
end

"""
     getRotationMatrix(axis::AbstractVector, angle::Real) --> SMatrix{3,3}

Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in radians.
Reference: [Rotation matrix from axis and angle](https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)

# Example

```julia
using LinearAlgebra
v = [-0.5, 1.0, 1.0]
v̂ = normalize(v)
θ = deg2rad(-74)
R = getRotationMatrix(v̂, θ)
```
"""
function getRotationMatrix(v::AbstractVector{<:AbstractFloat}, θ::Real)
   sinθ, cosθ = sincos(eltype(v)(θ))
   tmp = 1 - cosθ
   m = @SMatrix [cosθ+v[1]^2 * tmp v[1] * v[2] * tmp-v[3] * sinθ v[1] * v[3] * tmp+v[2] * sinθ;
                 v[1] * v[2] * tmp+v[3] * sinθ cosθ+v[2]^2 * tmp v[2] * v[3] * tmp-v[1] * sinθ;
                 v[1] * v[3] * tmp-v[2] * sinθ v[3] * v[2] * tmp+v[1] * sinθ cosθ+v[3]^2 * tmp]
end

"""
     generate_mock_amrex_data(output_dir::String; num_particles::Int=10, 
        real_component_names::Vector{String}=["u", "v"],
        particle_gen::Function)

Generate mock AMReX particle data for testing and benchmarking.

# Arguments

  - `output_dir::String`: Directory to save the mock data.
  - `num_particles::Int`: Number of particles to generate.
  - `real_component_names::Vector{String}`: Names of the extra real components (beyond x, y, z).
  - `particle_gen::Function`: A function `(i, n_reals) -> tuple` that takes an index `i` (1-based) and the total number of real components `n_reals`. It should return a tuple of `n_reals` Float64 values: `(x, y, z, comp1, comp2, ...)`.
"""
function generate_mock_amrex_data(output_dir::String;
      num_particles::Int = 10,
      real_component_names::Vector{String} = ["u", "v"],
      particle_gen::Function = (i, n_reals) -> (
         Float64(i), Float64(i), Float64(i), Float64(i * 10), Float64(i * 100))
)
   ptype = "particles"
   base_dir = joinpath(output_dir, ptype)
   mkpath(base_dir)

   n_extra = length(real_component_names)

   # Create Header
   header_path = joinpath(base_dir, "Header")
   open(header_path, "w") do f
      println(f, "Version_double")
      println(f, "3") # dim
      println(f, "$n_extra") # num_real_extra
      for name in real_component_names
         println(f, name)
      end
      println(f, "2") # num_int_extra (total int = 2 + 2 = 4)
      println(f, "id_1")
      println(f, "id_2")
      println(f, "0") # is_checkpoint (False)
      println(f, "$num_particles") # num_particles
      println(f, "$(num_particles + 1)") # max_next_id
      println(f, "0") # finest_level
      println(f, "1") # grids_per_level[0]
      println(f, "1 $num_particles 0") # grid info: which, count, where
   end

   # Create Level directory
   level_dir = joinpath(base_dir, "Level_0")
   mkpath(level_dir)

   # Create Particle_H
   particle_h_path = joinpath(level_dir, "Particle_H")
   open(particle_h_path, "w") do f
      println(f, "(1 0") # num_boxes level
      println(f, "((0,0,0) (10,10,10) (0,0,0))")
   end

   # Create Data file
   data_fn = joinpath(level_dir, "DATA_00001")
   open(data_fn, "w") do f
      # Write particles
      # structure: num_particles * (3+n_extra) reals
      # x, y, z, ...
      n_reals = 3 + n_extra
      data = zeros(Float64, num_particles * n_reals)
      for i in 1:num_particles
         vals = particle_gen(i, n_reals)
         for j in 1:n_reals
            data[(i - 1) * n_reals + j] = vals[j]
         end
      end
      write(f, data)
   end

   # Create Main Header (for domain info)
   main_header_path = joinpath(output_dir, "Header")
   open(main_header_path, "w") do f
      println(f, "HyperCLaw-V1.1")
      println(f, "0") # num_fields
      println(f, "3") # dim
      println(f, "0.0") # time
      println(f, "0") # refine_ratio
      println(f, "0.0 0.0 0.0") # left_edge
      println(f, "10.0 10.0 10.0") # right_edge
      println(f, "0")
      println(f, "((0,0,0) (10,10,10) (0,0,0))") # domain size
   end
end

function _triangulate_matplotlib end
