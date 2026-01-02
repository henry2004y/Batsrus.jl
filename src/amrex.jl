# AMReX particle data reader and analyzer.

using FHist

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
   ranges = (x_range, y_range, z_range)

   # Ensure data is loaded into memory. This will be a no-op if already loaded.
   # If real data is already loaded, filter in memory
   if !isnothing(data._rdata)
      rdata = data._rdata
      if isempty(rdata)
         return similar(rdata, size(rdata, 1), 0)
      end

      return _select_particles_in_memory(rdata, ranges, data.dim)
   end

   # If not loaded, optimize by reading specific grids
   return _select_particles_from_files(data, ranges)
end

function _select_particles_in_memory(rdata::Matrix{T}, ranges, dim) where T
   check_x = dim >= 1 && !isnothing(ranges[1])
   xlo, xhi = check_x ? ranges[1] : (T(0), T(0))

   check_y = dim >= 2 && !isnothing(ranges[2])
   ylo, yhi = check_y ? ranges[2] : (T(0), T(0))

   check_z = dim >= 3 && !isnothing(ranges[3])
   zlo, zhi = check_z ? ranges[3] : (T(0), T(0))

   valid_indices = filter(1:size(rdata, 2)) do i
      if check_x
         val = rdata[1, i]
         (val < xlo || val > xhi) && return false
      end
      if check_y
         val = rdata[2, i]
         (val < ylo || val > yhi) && return false
      end
      if check_z
         val = rdata[3, i]
         (val < zlo || val > zhi) && return false
      end
      return true
   end

   return rdata[:, valid_indices]
end

function _select_particles_from_files(data::AMReXParticle{T}, ranges) where T
   # Convert physical range to index range
   dx = SVector{3, Float64}(
      (data.right_edge[1] - data.left_edge[1]) / data.domain_dimensions[1],
      (data.right_edge[2] - data.left_edge[2]) / data.domain_dimensions[2],
      data.dim == 3 ? (data.right_edge[3] - data.left_edge[3]) / data.domain_dimensions[3] :
      1.0
   )

   target_idx_ranges = Vector{Union{Tuple{Int, Int}, Nothing}}(undef, data.dim)
   for i in 1:(data.dim)
      if !isnothing(ranges[i])
         idx_min = floor(Int, (ranges[i][1] - data.left_edge[i]) / dx[i])
         idx_max = floor(Int, (ranges[i][2] - data.left_edge[i]) / dx[i])
         target_idx_ranges[i] = (idx_min, idx_max)
      else
         target_idx_ranges[i] = nothing
      end
   end

   overlapping_grids = Tuple{Int, Int}[] # (level_num_0_indexed, grid_index_1_indexed)

   for (lvl_idx, boxes) in enumerate(data.level_boxes)
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

   selected_rdata = Vector{Matrix{T}}()

   base_fn = joinpath(data.output_dir, data.ptype)
   n_real = data.header.num_real
   dim = data.dim

   check_x = dim >= 1 && !isnothing(ranges[1])
   xlo, xhi = check_x ? ranges[1] : (T(0), T(0))

   check_y = dim >= 2 && !isnothing(ranges[2])
   ylo, yhi = check_y ? ranges[2] : (T(0), T(0))

   check_z = dim >= 3 && !isnothing(ranges[3])
   zlo, zhi = check_z ? ranges[3] : (T(0), T(0))

   for (level_num, grid_idx) in overlapping_grids
      # header.grids stores (which, count, where)
      grid_data = data.header.grids[level_num + 1][grid_idx]
      which, count, offset = grid_data

      if count == 0
         continue
      end

      data_fn = joinpath(base_fn, "Level_$(level_num)", @sprintf("DATA_%05d", which))

      open(data_fn, "r") do f
         seek(f, offset)
         if data.header.is_checkpoint
            # Skip integers
            skip(f, count * data.header.num_int * sizeof(data.header.int_type))
         end

         # Read floats
         floats_vec = Vector{T}(undef, count * n_real)
         read!(f, floats_vec)

         valid_rows = filter(0:(count - 1)) do k
            if check_x
               val = floats_vec[k * n_real + 1]
               (val < xlo || val > xhi) && return false
            end
            if check_y
               val = floats_vec[k * n_real + 2]
               (val < ylo || val > yhi) && return false
            end
            if check_z
               val = floats_vec[k * n_real + 3]
               (val < zlo || val > zhi) && return false
            end
            return true
         end

         if !isempty(valid_rows)
            # Create result matrix for this block
            res = Matrix{T}(undef, n_real, length(valid_rows))
            for (i, row_idx) in enumerate(valid_rows)
               base = row_idx * n_real
               for j in 1:n_real
                  res[j, i] = floats_vec[base + j]
               end
            end
            push!(selected_rdata, res)
         end
      end
   end

   if isempty(selected_rdata)
      return Matrix{T}(undef, n_real, 0)
   end

   return hcat(selected_rdata...)
end

const _ALIAS_MAP = Dict(
   "vx" => "velocity_x",
   "vy" => "velocity_y",
   "vz" => "velocity_z"
)

_resolve_alias(variable_name::String) = get(_ALIAS_MAP, variable_name, variable_name)

"""
    get_phase_space_density(data, x_var, y_var; bins=100, x_range=nothing, y_range=nothing, z_range=nothing)::Hist2D

Calculates the 2D phase space density for selected variables.
"""
function get_phase_space_density(
      data::AMReXParticle{T},
      x_variable::String,
      y_variable::String;
      bins::Union{Int, Tuple{Int, Int}} = 100,
      x_range = nothing,
      y_range = nothing,
      z_range = nothing,
      transform::Union{Function, Nothing} = nothing,
      normalize::Bool = false
) where T
   # Select data
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
      rdata = new_data

      # Update component mapping for variable lookup below
      component_names = new_names
      component_map = Dict(name => i for (i, name) in enumerate(component_names))
   else
      # Map component names to columns (original)
      component_names = data.header.real_component_names
      component_map = Dict(name => i for (i, name) in enumerate(component_names))
   end

   x_variable = _resolve_alias(x_variable)
   y_variable = _resolve_alias(y_variable)

   if !haskey(component_map, x_variable) || !haskey(component_map, y_variable)
      error("Invalid variable name. Available: $(keys(component_map))")
   end

   x_index = component_map[x_variable]
   y_index = component_map[y_variable]

   x_data = rdata[x_index, :]
   y_data = rdata[y_index, :] # FHist uses nbins keyword
   arg_bins = bins isa Int ? (bins, bins) : bins
   nx, ny = arg_bins

   # Calculate edges explicitly
   if !isnothing(x_range)
      xmin, xmax = x_range
   else
      xmin, xmax = extrema(x_data)
   end

   if !isnothing(y_range)
      ymin, ymax = y_range
   else
      ymin, ymax = extrema(y_data)
   end

   x_edges = range(xmin, xmax, length = nx + 1)
   y_edges = range(ymin, ymax, length = ny + 1)

   h = Hist2D((x_data, y_data); binedges = (x_edges, y_edges))

   # Normalize to probability density if requested
   if normalize
      h = FHist.normalize(h)
   end

   return h
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
   # Try explicit names first, then aliases
   possible_names = [
      ["vx", "vy", "vz"], ["ux", "uy", "uz"], ["velocity_x", "velocity_y", "velocity_z"]]
   vel_indices = Int[]

   component_names = data.header.real_component_names
   component_map = Dict(name => i for (i, name) in enumerate(component_names))

   # Find valid velocity columns
   for names in possible_names
      indices = [get(component_map, n, get(component_map, _resolve_alias(n), 0))
                 for n in names[1:vdim]]
      if all(i -> i > 0, indices)
         vel_indices = indices
         break
      end
   end

   if isempty(vel_indices)
      error("Could not identify velocity components for vdim=$vdim. Checked standard names (v, u, velocity).")
   end

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
