# AMReX particle data reader and analyzer.

using FHist
using GaussianMixtures

struct AMReXParticleHeader
   version_string::String
   real_type::DataType
   int_type::DataType
   dim::Int
   num_int_base::Int
   num_real_base::Int
   real_component_names::Vector{String}
   int_component_names::Vector{String}
   num_real_extra::Int
   num_int_extra::Int
   num_int::Int
   num_real::Int
   is_checkpoint::Bool
   num_particles::Int
   max_next_id::Int
   finest_level::Int
   num_levels::Int
   grids_per_level::Vector{Int}
   # grids[level][grid_index] = (which, count, offset)
   grids::Vector{Vector{Tuple{Int, Int, Int}}}
end

function AMReXParticleHeader(header_filename::AbstractString)
   grids = Vector{Vector{Tuple{Int, Int, Int}}}()
   grids_per_level = Int[]
   int_component_names = String[]
   real_component_names = String[]

   open(header_filename, "r") do f
      version_string = readline(f) |> strip
      particle_real_type = split(version_string, '_')[end]
      if particle_real_type == "double"
         real_type = Float64
      elseif particle_real_type == "single"
         real_type = Float32
      else
         error("Did not recognize particle real type: $particle_real_type")
      end
      int_type = Int32

      dim = parse(Int, readline(f))
      num_int_base = 2
      num_real_base = dim

      if dim == 3
         real_component_names = ["x", "y", "z"]
      elseif dim == 2
         real_component_names = ["x", "y"]
      end

      int_component_names = ["particle_id", "particle_cpu"]

      num_real_extra = parse(Int, readline(f))
      for _ in 1:num_real_extra
         push!(real_component_names, readline(f) |> strip)
      end

      num_int_extra = parse(Int, readline(f))
      for _ in 1:num_int_extra
         push!(int_component_names, readline(f) |> strip)
      end

      num_int = num_int_base + num_int_extra
      num_real = num_real_base + num_real_extra

      is_checkpoint = parse(Int, readline(f)) |> Bool
      num_particles = parse(Int, readline(f))
      max_next_id = parse(Int, readline(f))
      finest_level = parse(Int, readline(f))
      num_levels = finest_level + 1

      if !is_checkpoint
         num_int_base = 0
         num_int_extra = 0
         num_int = 0
      end

      grids_per_level = zeros(Int, num_levels)
      for level_num in 1:num_levels
         grids_per_level[level_num] = parse(Int, readline(f))
      end

      grids = [Vector{Tuple{Int, Int, Int}}() for _ in 1:num_levels]

      for level_num in 1:num_levels
         for _ in 1:grids_per_level[level_num]
            entry = [parse(Int, val) for val in split(readline(f))]
            push!(grids[level_num], Tuple(entry))
         end
      end

      AMReXParticleHeader(
         version_string,
         real_type,
         int_type,
         dim,
         num_int_base,
         num_real_base,
         real_component_names,
         int_component_names,
         num_real_extra,
         num_int_extra,
         num_int,
         num_real,
         is_checkpoint,
         num_particles,
         max_next_id,
         finest_level,
         num_levels,
         grids_per_level,
         grids
      )
   end
end

function Base.show(io::IO, header::AMReXParticleHeader)
   println(io, "Version string: ", header.version_string)
   println(io, "Dimensions: ", header.dim)
   println(io, "Number of integer components: ", header.num_int)
   println(io, "Integer component names: ", header.int_component_names)
   println(io, "Number of real components: ", header.num_real)
   println(io, "Real component names: ", header.real_component_names)
   println(io, "Is checkpoint: ", header.is_checkpoint)
   println(io, "Number of particles: ", header.num_particles)
   println(io, "Max next ID: ", header.max_next_id)
   println(io, "Finest level: ", header.finest_level)
   println(io, "Number of levels: ", header.num_levels)
   for level_num in 1:(header.num_levels)
      println(
         io, "  Level ", level_num - 1, ": ", header.grids_per_level[level_num], " grids")
   end
end

mutable struct AMReXParticle{T <: Real}
   const output_dir::String
   const ptype::String
   _idata::Union{Matrix{Int32}, Nothing}
   _rdata::Union{Matrix{T}, Nothing}
   const level_boxes::Vector{Vector{Tuple{
      Tuple{Int, Vararg{Int}}, Tuple{Int, Vararg{Int}}}}}
   const header::AMReXParticleHeader
   const dim::Int
   const time::Float64
   const left_edge::Vector{Float64}
   const right_edge::Vector{Float64}
   const domain_dimensions::Vector{Int}

   function AMReXParticle(output_dir::AbstractString)
      ptype = "particles"
      idata = nothing
      rdata = nothing

      # Parse main header
      header_path = joinpath(output_dir, "Header")

      dim, time, left_edge, right_edge, domain_dimensions = open(header_path, "r") do f
         readline(f) # version string
         num_fields = parse(Int, readline(f))
         for _ in 1:num_fields
            readline(f)
         end

         dim_val = parse(Int, readline(f))
         time_val = parse(Float64, readline(f))
         readline(f) # prob_refine_ratio

         left_edge_val = [parse(Float64, v) for v in split(readline(f))]
         right_edge_val = [parse(Float64, v) for v in split(readline(f))]
         readline(f)

         dim_line = readline(f) |> strip
         matches = [parse(Int, m.match) for m in eachmatch(r"-?\d+", dim_line)]

         local domain_dimensions_val
         if length(matches) >= 2 * dim_val
            coords = matches[1:(2 * dim_val)]
            domain_dimensions_val = [coords[i + dim_val] - coords[i] + 1 for i in 1:dim_val]
         else
            error("Dimension mismatch when parsing domain dimensions from Header.")
         end

         (dim_val, time_val, left_edge_val, right_edge_val, domain_dimensions_val)
      end

      header = AMReXParticleHeader(joinpath(output_dir, ptype, "Header"))
      level_boxes = _read_level_boxes(output_dir, ptype, header)

      T = header.real_type

      new{T}(output_dir, ptype, idata, rdata, level_boxes, header,
         dim, time, left_edge, right_edge, domain_dimensions)
   end
end

function _read_level_boxes(output_dir::String, ptype::String, header::AMReXParticleHeader)
   level_boxes = [Vector{Tuple{Tuple{Int, Vararg{Int}}, Tuple{Int, Vararg{Int}}}}()
                  for _ in 1:(header.num_levels)]

   for level_num in 0:(header.num_levels - 1)
      particle_h_path = joinpath(
         output_dir, ptype, "Level_$(level_num)", "Particle_H")

      if !isfile(particle_h_path)
         continue
      end

      lines = readlines(particle_h_path)

      boxes = Vector{Tuple{Tuple{Int, Vararg{Int}}, Tuple{Int, Vararg{Int}}}}()

      @views for line in lines[2:end] # Skip first line
         line = strip(line)
         if startswith(line, "((") && endswith(line, "))")
            parts = [parse(Int, m.match) for m in eachmatch(r"-?\d+", line)]

            if header.dim == 2 && length(parts) >= 4
               lo_corner = (parts[1], parts[2])
               hi_corner = (parts[3], parts[4])
               push!(boxes, (lo_corner, hi_corner))
            elseif header.dim == 3 && length(parts) >= 6
               lo_corner = (parts[1], parts[2], parts[3])
               hi_corner = (parts[4], parts[5], parts[6])
               push!(boxes, (lo_corner, hi_corner))
            end
         end
      end
      level_boxes[level_num + 1] = boxes
   end

   return level_boxes
end

function read_amrex_binary_particle_file(fn::AbstractString, header::AMReXParticleHeader)
   ptype = "particles"
   base_fn = joinpath(fn, ptype)

   idata = Matrix{header.int_type}(undef, header.num_int, header.num_particles)
   rdata = Matrix{header.real_type}(undef, header.num_real, header.num_particles)

   ip = 1
   for lvl in 0:(header.num_levels - 1)
      level_grids = header.grids[lvl + 1]
      for (which, count, offset) in level_grids
         if count == 0
            continue
         end

         data_fn = joinpath(base_fn, "Level_$(lvl)", @sprintf("DATA_%05d", which))

         open(data_fn, "r") do f
            seek(f, offset)
            if header.is_checkpoint
               # Read integers
               ints_vec = Vector{header.int_type}(undef, count * header.num_int)
               read!(f, ints_vec)
               # Reshape and assign.
               # AMReX stores particles contiguously (AoS). Reading linearly and reshaping
               # to (num_int, count) preserves the particle-major order in memory (SoA in Julia).
               ints_mat = reshape(ints_vec, header.num_int, count)
               idata[:, ip:(ip + count - 1)] = ints_mat
            end

            # Read floats
            floats_vec = Vector{header.real_type}(undef, count * header.num_real)
            read!(f, floats_vec)
            # Reshape to (num_real, count) where columns are particles.
            floats_mat = reshape(floats_vec, header.num_real, count)
            rdata[:, ip:(ip + count - 1)] = floats_mat
         end
         ip += count
      end
   end

   return idata, rdata
end

function load_data!(data::AMReXParticle{T}) where T
   if isnothing(data._idata) || isnothing(data._rdata)
      idata, rdata = read_amrex_binary_particle_file(data.output_dir, data.header)
      data._idata = idata
      data._rdata = rdata
   end
end

function get_idata(data::AMReXParticle{T}) where T
   load_data!(data)
   return data._idata
end

function get_rdata(data::AMReXParticle{T}) where T
   load_data!(data)
   return data._rdata
end

Base.getproperty(obj::AMReXParticle{T}, sym::Symbol) where T =
   if sym === :idata
      return get_idata(obj)
   elseif sym === :rdata
      return get_rdata(obj)
   else
      return getfield(obj, sym)
   end

function Base.show(io::IO, data::AMReXParticle{T}) where T
   println(io, "AMReXParticle{$T} from ", data.output_dir)
   println(io, "Time: ", data.time)
   println(io, "Dimensions: ", data.dim)
   println(io, "Domain Dimensions: ", data.domain_dimensions)
   println(io, "Domain Edges: ", data.left_edge, " to ", data.right_edge)
   println(io, "Integer component names: ", data.header.int_component_names)
   println(io, "Real component names: ", data.header.real_component_names)

   if !isnothing(data._idata)
      println(io, "Particle data shape (int): ", size(data._idata))
      println(io, "Particle data shape (real): ", size(data._rdata))
   end
end

function select_particles_in_region(
      data::AMReXParticle{T};
      x_range = nothing,
      y_range = nothing,
      z_range = nothing
) where T
   # Ensure data is loaded into memory.
   # If real data is already loaded, filter in memory
   rdata = data._rdata
   if !isnothing(rdata)
      if isempty(rdata)
         return similar(rdata, size(rdata, 1), 0)
      end

      return _select_particles_in_memory(rdata, x_range, y_range, z_range, data.dim)
   end

   # If not loaded, optimize by reading specific grids
   return _select_particles_from_files(data, x_range, y_range, z_range)
end

function _select_particles_in_memory(
      rdata::Matrix{T}, x_range, y_range, z_range, dim) where T
   check_x = dim >= 1 && !isnothing(x_range)
   xlo, xhi = check_x ? x_range : (T(0), T(0))

   check_y = dim >= 2 && !isnothing(y_range)
   ylo, yhi = check_y ? y_range : (T(0), T(0))

   check_z = dim >= 3 && !isnothing(z_range)
   zlo, zhi = check_z ? z_range : (T(0), T(0))

   n_vars, n_particles = size(rdata)

   # Pass 1: Count valid particles to avoid vector allocation
   count = 0
   @inbounds for i in 1:n_particles
      keep = true
      if check_x
         val = rdata[1, i]
         if val < xlo || val > xhi
            keep = false
         end
      end
      if keep && check_y
         val = rdata[2, i]
         if val < ylo || val > yhi
            keep = false
         end
      end
      if keep && check_z
         val = rdata[3, i]
         if val < zlo || val > zhi
            keep = false
         end
      end
      if keep
         count += 1
      end
   end

   # Allocate result matrix exactly
   new_data = Matrix{T}(undef, n_vars, count)

   # Pass 2: Fill data
   current_idx = 0
   @inbounds for i in 1:n_particles
      keep = true
      if check_x
         val = rdata[1, i]
         if val < xlo || val > xhi
            keep = false
         end
      end
      if keep && check_y
         val = rdata[2, i]
         if val < ylo || val > yhi
            keep = false
         end
      end
      if keep && check_z
         val = rdata[3, i]
         if val < zlo || val > zhi
            keep = false
         end
      end

      if keep
         current_idx += 1
         for k in 1:n_vars
            new_data[k, current_idx] = rdata[k, i]
         end
      end
   end

   return new_data
end

function _select_particles_from_files(
      data::AMReXParticle{T}, x_range, y_range, z_range) where T
   ranges = (x_range, y_range, z_range)
   # Convert physical range to index range
   dx = SVector{3, Float64}(
      (data.right_edge[1] - data.left_edge[1]) / data.domain_dimensions[1],
      (data.right_edge[2] - data.left_edge[2]) / data.domain_dimensions[2],
      data.dim == 3 ? (data.right_edge[3] - data.left_edge[3]) / data.domain_dimensions[3] :
      1.0
   )

   target_idx_ranges = Vector{Union{Tuple{Int, Int}, Nothing}}(undef, data.dim)
   @inbounds for i in 1:(data.dim)
      if !isnothing(ranges[i])
         idx_min = floor(Int, (ranges[i][1] - data.left_edge[i]) / dx[i])
         idx_max = floor(Int, (ranges[i][2] - data.left_edge[i]) / dx[i])
         target_idx_ranges[i] = (idx_min, idx_max)
      else
         target_idx_ranges[i] = nothing
      end
   end

   overlapping_grids = Tuple{Int, Int}[] # (level_num_0_indexed, grid_index_1_indexed)

   @inbounds for (lvl_idx, boxes) in enumerate(data.level_boxes)
      level_num = lvl_idx - 1
      for (grid_idx, (lo, hi)) in enumerate(boxes)
         box_overlap = true
         for i in 1:(data.dim)
            if !isnothing(target_idx_ranges[i])
               target_min, target_max = target_idx_ranges[i]
               # lo and hi are tuples of Int indices
               if hi[i] < target_min || lo[i] > target_max
                  box_overlap = false
                  break
               end
            end
         end
         if box_overlap
            push!(overlapping_grids, (level_num, grid_idx))
         end
      end
   end

   # Preallocate result container (max possible size = number of grids)
   # We use a counter to track actual filled grids
   selected_rdata = Vector{Matrix{T}}(undef, length(overlapping_grids))
   filled_count = 0

   n_real = data.header.num_real
   base_fn = joinpath(data.output_dir, data.ptype)

   # Function barrier to dispatch on n_real for SVector optimization
   filled_count = _process_grids!(
      selected_rdata,
      Val(n_real),
      overlapping_grids,
      data,
      base_fn,
      x_range,
      y_range,
      z_range
   )

   if filled_count == 0
      return Matrix{T}(undef, n_real, 0)
   end

   resize!(selected_rdata, filled_count)
   return reduce(hcat, selected_rdata)
end

function _process_grids!(
      selected_rdata::Vector{Matrix{T}},
      ::Val{N},
      overlapping_grids,
      data,
      base_fn,
      x_range,
      y_range,
      z_range
) where {T, N}
   filled_idx = 0

   dim = data.dim

   check_x = dim >= 1 && !isnothing(x_range)
   xlo, xhi = check_x ? x_range : (T(0), T(0))

   check_y = dim >= 2 && !isnothing(y_range)
   ylo, yhi = check_y ? y_range : (T(0), T(0))

   check_z = dim >= 3 && !isnothing(z_range)
   zlo, zhi = check_z ? z_range : (T(0), T(0))

   @inbounds for (level_num, grid_idx) in overlapping_grids
      # header.grids stores (which, count, where)
      grid_data = data.header.grids[level_num + 1][grid_idx]
      which, count, offset = grid_data

      count == 0 && continue

      data_fn = joinpath(base_fn, "Level_$(level_num)", @sprintf("DATA_%05d", which))

      open(data_fn, "r") do f
         seek(f, offset)
         if data.header.is_checkpoint
            # Skip integers
            skip(f, count * data.header.num_int * sizeof(data.header.int_type))
         end

         # Read floats
         floats_vec = Vector{T}(undef, count * N)
         read!(f, floats_vec)

         # Reinterpret as SVector for fast access
         vectors = reinterpret(SVector{N, T}, floats_vec)

         # Pass 1: Count valid
         valid_count = 0
         @inbounds for k in 1:count
            val = vectors[k]
            keep = true
            if check_x
               v = val[1]
               if v < xlo || v > xhi
                  keep = false
               end
            end
            if keep && check_y
               v = val[2]
               if v < ylo || v > yhi
                  keep = false
               end
            end
            if keep && check_z
               v = val[3]
               if v < zlo || v > zhi
                  keep = false
               end
            end
            if keep
               valid_count += 1
            end
         end

         if valid_count > 0
            # Pass 2: Fill
            # Allocate result matrix for this grid
            res = Matrix{T}(undef, N, valid_count)

            curr = 0
            @inbounds for k in 1:count
               val = vectors[k]
               keep = true
               if check_x
                  v = val[1]
                  if v < xlo || v > xhi
                     keep = false
                  end
               end
               if keep && check_y
                  v = val[2]
                  if v < ylo || v > yhi
                     keep = false
                  end
               end
               if keep && check_z
                  v = val[3]
                  if v < zlo || v > zhi
                     keep = false
                  end
               end

               if keep
                  curr += 1
                  res[:, curr] = val
               end
            end

            filled_idx += 1
            selected_rdata[filled_idx] = res
         end
      end
   end

   return filled_idx
end

const _ALIAS_MAP = Dict(
   "vx" => "velocity_x",
   "vy" => "velocity_y",
   "vz" => "velocity_z"
)

_resolve_alias(variable_name::String) = get(_ALIAS_MAP, variable_name, variable_name)

"""
    get_phase_space_density(data, variables...; bins=100, x_range=nothing, y_range=nothing, z_range=nothing, edges=nothing)::Hist

Calculates the phase space density for selected variables.
Supports 1D, 2D, and 3D histograms.

# Arguments

  - `data`: AMReXParticle data object.
  - `variables`: Variable names to compute the histogram for (e.g., "vx", "vy").
  - `bins`: Number of bins for the histogram.
  - `x_range`, `y_range`, `z_range`: **Spatial selection ranges**. Only particles within these ranges in configuration space are selected.
  - `edges`: **Histogram binning edges**. If provided, these define the exact bins for the `variables`. If not provided, bins are determined automatically from the data extrema.
  - `transform`: Optional function to transform the data before binning.
  - `normalize`: Whether to normalize the histogram to a probability density (default: `false`).
"""
function get_phase_space_density(
      data::AMReXParticle{T},
      variables::Vararg{String, N};
      bins::Union{Int, Tuple{Vararg{Int, N}}} = 100,
      edges = nothing,
      x_range = nothing,
      y_range = nothing,
      z_range = nothing,
      transform::Union{Function, Nothing} = nothing,
      normalize::Bool = false
) where {T, N}
   # Select data
   local rdata::Matrix{T}
   if !isnothing(x_range) || !isnothing(y_range) || !isnothing(z_range)
      rdata = select_particles_in_region(data; x_range, y_range, z_range)
   else
      rdata = data.rdata
   end

   if isempty(rdata)
      error("No particles found for phase space density calculation.")
   end

   if !isnothing(transform)
      # Apply transform
      # Need full names including x, y, z
      names = data.header.real_component_names
      full_names = ["x", "y", "z", names...]

      # rdata from select_particles_in_region has ALL real components.
      new_data, new_names = transform(rdata, full_names)

      # Use transformed data
      rdata = new_data::Matrix{T}

      # Update component mapping for variable lookup below
      component_names = new_names::Vector{String}
      component_map = Dict(name => i for (i, name) in enumerate(component_names))
   else
      # Map component names to columns (original)
      component_names = data.header.real_component_names
      component_map = Dict(name => i for (i, name) in enumerate(component_names))
   end

   for var in variables
      var = _resolve_alias(var)
      if !haskey(component_map, var)
         error("Invalid variable name: $var. Available: $(keys(component_map))")
      end
   end

   # Check for weights
   weight_index = get(component_map, "weight", 0)
   weights = weight_index > 0 ? rdata[weight_index, :] : nothing

   indices = ntuple(i -> component_map[_resolve_alias(variables[i])], N)
   selected_data = ntuple(i -> rdata[indices[i], :], N)

   # Handle bins
   if bins isa Int
      arg_bins = ntuple(_ -> bins, N)
   else
      arg_bins = bins
      if length(arg_bins) != N
         error("Length of bins must match number of variables")
      end
   end

   edges = ntuple(
      i -> begin
         if !isnothing(edges) && i <= length(edges) && !isnothing(edges[i])
            return edges[i]
         end

         data_i = selected_data[i]
         nbins = arg_bins[i]

         # Use explicit limits if provided, otherwise data extrema
         vmin, vmax = extrema(data_i)

         # Handle case of single value or empty
         if vmin == vmax
            vmin -= 0.5
            vmax += 0.5
         else
            # Ensure the max value is included in the last bin
            vmax = nextfloat(vmax)
         end

         range(vmin, vmax, length = nbins + 1)
      end, N)

   h = _create_hist(selected_data, edges, weights)

   return _finalize_hist(h, Val(normalize))
end

function _create_hist(selected_data, edges, weights)
   N = length(selected_data)
   if N == 1
      Hist1D(selected_data[1]; binedges = edges[1], weights)
   elseif N == 2
      Hist2D(Tuple(selected_data); binedges = edges, weights)
   elseif N == 3
      Hist3D(Tuple(selected_data); binedges = edges, weights)
   else
      error("Only 1D, 2D, and 3D phase space densities are supported.")
   end
end

_finalize_hist(h, ::Val{true}) = FHist.normalize(h)

_finalize_hist(h, ::Val{false}) = h

function _get_velocity_indices(data::AMReXParticle{T}, vdim::Int) where T
   possible_names = (
      ("vx", "vy", "vz"), ("ux", "uy", "uz"), ("velocity_x", "velocity_y", "velocity_z"))

   component_names = data.header.real_component_names
   component_map = Dict{String, Int}(name => i for (i, name) in enumerate(component_names))

   vel_indices = Vector{Int}(undef, vdim)

   for names in possible_names
      found_all = true
      for k in 1:vdim
         nm = names[k]
         idx = get(component_map, nm, get(component_map, _resolve_alias(nm), 0))
         if idx == 0
            found_all = false
            break
         end
         vel_indices[k] = idx
      end

      if found_all
         return vel_indices
      end
   end

   error("Could not identify velocity components for vdim=$vdim. Checked standard names (v, u, velocity).")
end

"""
    classify_particles(data, region; vdim=3, bulk_vel=nothing, vth=nothing, nsigma=3.0)

Classify particles in a spatial region into core Maxwellian and suprathermal populations.

# Arguments

  - `data::AMReXParticle`: Particle data.
  - `region`: Passed as kwargs `x_range`, `y_range`, `z_range`.
  - `vdim`: Velocity dimension (1, 2, or 3).
  - `bulk_vel`: Core bulk velocity. If `nothing`, estimated from peak density.
  - `vth`: Core thermal velocity. Must be provided.
  - `nsigma`: Threshold for classification in units of thermal velocity.

# Returns

  - `(core, suprathermal)`: Two matrices containing the classified particles.
"""
function classify_particles(
      data::AMReXParticle{T};
      x_range = nothing,
      y_range = nothing,
      z_range = nothing,
      vdim::Int = 3,
      bulk_vel = nothing,
      vth = nothing,
      nsigma = 3.0
) where T
   # 1. Select particles
   particles = select_particles_in_region(data; x_range, y_range, z_range)
   if isempty(particles)
      return particles, particles
   end

   # 2. Identify velocity columns
   vel_indices = _get_velocity_indices(data, vdim)

   velocities = particles[vel_indices, :]

   # 3. Determine bulk velocity
   if isnothing(bulk_vel)
      nbins = 50
      detected_bulk = zeros(T, vdim)

      if vdim == 1
         h = Hist1D(velocities[1, :], nbins = nbins)
         _, max_idx = findmax(h.bincounts)
         edges = h.binedges isa Tuple ? h.binedges[1] : h.binedges
         detected_bulk[1] = (edges[max_idx] + edges[max_idx + 1]) / 2
      elseif vdim == 2
         h = Hist2D((velocities[1, :], velocities[2, :]), nbins = (nbins, nbins))
         _, max_idx = findmax(h.bincounts)
         # max_idx is CartesianIndex
         x_edges = h.binedges[1]
         y_edges = h.binedges[2]
         detected_bulk[1] = (x_edges[max_idx[1]] + x_edges[max_idx[1] + 1]) / 2
         detected_bulk[2] = (y_edges[max_idx[2]] + y_edges[max_idx[2] + 1]) / 2
      elseif vdim == 3
         # Use marginal peaks for 3D as a robust fallback
         for i in 1:vdim
            h = Hist1D(velocities[i, :], nbins = nbins)
            _, max_idx = findmax(h.bincounts)
            edges = h.binedges isa Tuple ? h.binedges[1] : h.binedges
            detected_bulk[i] = (edges[max_idx] + edges[max_idx + 1]) / 2
         end
      end
      bulk_vel = detected_bulk
   end

   # 4. Thermal velocity handling
   if isnothing(vth)
      error("Thermal velocity `vth` must be provided.")
   end

   if vth isa Real
      vth_vec = fill(T(vth), vdim)
   else
      vth_vec = vth
   end

   # 5. Classify
   n_part = size(particles, 2)
   is_core = Vector{Bool}(undef, n_part)
   threshold_sq = nsigma^2

   for i in 1:n_part
      d2 = zero(T)
      for k in 1:vdim
         d2 += ((velocities[k, i] - bulk_vel[k]) / vth_vec[k])^2
      end
      is_core[i] = d2 <= threshold_sq
   end

   return particles[:, is_core], particles[:, .!is_core]
end

"""
    fit_particle_velocity_gmm(data, n_clusters; x_range=nothing, y_range=nothing, z_range=nothing, vdim=3)

Fit a Gaussian Mixture Model to particle velocities in a region.

# Arguments

  - `data`: AMReXParticle data.
  - `n_clusters`: Number of GMM components.
  - `vdim`: Velocity dimension (1, 2, or 3).

# Returns

  - A vector of named tuples sorted by weight, each containing:

      + `weight`: Component weight.
      + `mean`: Component mean velocity (vector of length vdim).
      + `vth`: Component thermal velocity (vector of length vdim).
"""
function fit_particle_velocity_gmm(
      data::AMReXParticle{T},
      n_clusters::Int;
      x_range = nothing,
      y_range = nothing,
      z_range = nothing,
      vdim::Int = 3
) where T
   # 1. Select particles
   particles = select_particles_in_region(data; x_range, y_range, z_range)
   if isempty(particles)
      error("No particles found in the specified region.")
   end

   # 2. Extract velocities
   vel_indices = _get_velocity_indices(data, vdim)

   velocities = particles[vel_indices, :] # (vdim, n_particles)

   # 3. Fit GMM
   # GaussianMixtures expects (n_samples, n_features). It supports Float32/Float64.
   X = Matrix{T}(velocities')

   # kind=:diag for diagonal covariance (independent velocity components)
   gmm = GMM(n_clusters, X, kind = :diag)

   # 4. Interpret results
   results = [(
                 weight = T(gmm.w[i]),
                 mean = T.(gmm.μ[i, :]),
                 vth = T.(sqrt.(2 .* gmm.Σ[i, :]))
              ) for i in 1:n_clusters]

   # Sort by weight descending
   sort!(results, by = x -> x.weight, rev = true)

   return results
end
