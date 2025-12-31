# AMReX particle data reader and analyzer.

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
   # grids[level][grid_index] = (which, count, where)
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

mutable struct AMReXParticleData
   const output_dir::String
   const ptype::String
   _idata::Union{Matrix{Int32}, Nothing}
   _rdata::Union{Matrix{Float64}, Matrix{Float32}, Nothing}
   const level_boxes::Vector{Vector{Tuple{Tuple{Int, Vararg{Int}}, Tuple{Int, Vararg{Int}}}}}
   const header::AMReXParticleHeader
   const dim::Int
   const time::Float64
   const left_edge::Vector{Float64}
   const right_edge::Vector{Float64}
   const domain_dimensions::Vector{Int}

   function AMReXParticleData(output_dir::AbstractString)
      ptype = "particles"
      idata = nothing
      rdata = nothing
      
      # Parse main header
      header_path = joinpath(output_dir, "Header")

      local dim, time, left_edge, right_edge, domain_dimensions

      open(header_path, "r") do f
         readline(f) # version string
         num_fields = parse(Int, readline(f))
         for _ in 1:num_fields
            readline(f)
         end

         dim = parse(Int, readline(f))
         time = parse(Float64, readline(f))
         readline(f) # prob_refine_ratio

         left_edge = [parse(Float64, v) for v in split(readline(f))]
         right_edge = [parse(Float64, v) for v in split(readline(f))]
         readline(f)

         dim_line = readline(f) |> strip
         matches = [parse(Int, m.match) for m in eachmatch(r"\d+", dim_line)]

         # Let's adapt to dynamic dimensions
         if length(matches) >= 2 * dim
            coords = matches[1:(2*dim)]
            domain_dimensions = [coords[i+dim] - coords[i] + 1 for i in 1:dim]
         else
            error("Dimension mismatch when parsing domain dimensions from Header.")
         end
      end

      header = AMReXParticleHeader(joinpath(output_dir, ptype, "Header"))
      level_boxes = _read_level_boxes(output_dir, ptype, header)

      new(output_dir, ptype, idata, rdata, level_boxes, header,
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

      for line in lines[2:end] # Skip first line
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

   idata = Matrix{header.int_type}(undef, header.num_particles, header.num_int)
   rdata = Matrix{header.real_type}(undef, header.num_particles, header.num_real)

   ip = 1
   for lvl in 0:(header.num_levels - 1)
      level_grids = header.grids[lvl + 1]
      for (which, count, where) in level_grids
         if count == 0
            continue
         end

         data_fn = joinpath(base_fn, "Level_$(lvl)", @sprintf("DATA_%05d", which))

         open(data_fn, "r") do f
            seek(f, where)
            if header.is_checkpoint
               # Read integers
               ints_vec = Vector{header.int_type}(undef, count * header.num_int)
               read!(f, ints_vec)
               # Reshape and assign. Note: Julia matrices are column-major, but file is likely row-major (C-style).
               # AMReX stores particles contiguously, structure of array (N, num_int).
               # (num_int, count) -> (count, num_int)
               ints_mat = reshape(ints_vec, header.num_int, count)'
               idata[ip:(ip + count - 1), :] = ints_mat
            end

            # Read floats
            floats_vec = Vector{header.real_type}(undef, count * header.num_real)
            read!(f, floats_vec)
            # Reshape assuming C-order row-major storage for particles
            floats_mat = reshape(floats_vec, header.num_real, count)'
            rdata[ip:(ip + count - 1), :] = floats_mat
         end
         ip += count
      end
   end

   return idata, rdata
end

function load_data!(data::AMReXParticleData)
   if isnothing(data._idata) || isnothing(data._rdata)
      idata, rdata = read_amrex_binary_particle_file(data.output_dir, data.header)
      data._idata = idata
      data._rdata = rdata
   end
end

function get_idata(data::AMReXParticleData)
   load_data!(data)
   return data._idata
end

function get_rdata(data::AMReXParticleData)
   load_data!(data)
   return data._rdata
end

Base.getproperty(obj::AMReXParticleData, sym::Symbol) =
   if sym === :idata
      return get_idata(obj)
   elseif sym === :rdata
      return get_rdata(obj)
   else
      return getfield(obj, sym)
   end

function Base.show(io::IO, data::AMReXParticleData)
   println(io, "AMReXParticleData from ", data.output_dir)
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
      data::AMReXParticleData;
      x_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      y_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      z_range::Union{Tuple{Float64, Float64}, Nothing} = nothing
)
   # Ensure data is loaded into memory. This will be a no-op if already loaded.
   load_data!(data)
   rdata = data._rdata

   if isnothing(rdata)
      return Matrix{data.header.real_type}(undef, 0, data.header.num_real)
   end

   if isempty(rdata)
      return similar(rdata, 0, size(rdata, 2))
   end

   mask = trues(size(rdata, 1))
   ranges = [x_range, y_range, z_range]

   for i in 1:min(data.dim, length(ranges))
      if !isnothing(ranges[i])
         # AMReX real components are x, y, z, ...
         # so column i corresponds to the i-th dimension.
         col_data = @view rdata[:, i]
         mask .&= (col_data .>= ranges[i][1]) .& (col_data .<= ranges[i][2])
      end
   end

   return rdata[mask, :]
end
