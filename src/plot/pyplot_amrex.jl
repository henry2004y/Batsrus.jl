# Plotting functionalities for AMReX particle data.

using FHist

export get_phase_space_density

const _AXIS_LABEL_MAP = Dict(
   "velocity_x" => "\$v_x\$",
   "velocity_y" => "\$v_y\$",
   "velocity_z" => "\$v_z\$"
)

const _ALIAS_MAP = Dict(
   "vx" => "velocity_x",
   "vy" => "velocity_y",
   "vz" => "velocity_z"
)

_resolve_alias(variable_name::String) = get(_ALIAS_MAP, variable_name, variable_name)

_get_axis_label(variable_name::String) = get(_AXIS_LABEL_MAP, variable_name, variable_name)

"""
    get_phase_space_density(data, x_var, y_var; bins=100, x_range=nothing, y_range=nothing, z_range=nothing)::(H, xedges, yedges)

Calculates the 2D phase space density for selected variables.
"""
function get_phase_space_density(
      data::AMReXParticleData,
      x_variable::String,
      y_variable::String;
      bins::Union{Int, Tuple{Int, Int}} = 100,
      x_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      y_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      z_range::Union{Tuple{Float64, Float64}, Nothing} = nothing
)
   # Select data
   if !isnothing(x_range) || !isnothing(y_range) || !isnothing(z_range)
      rdata = select_particles_in_region(data; x_range, y_range, z_range)
   else
      rdata = data.rdata
   end

   if isempty(rdata)
      error("No particles found for phase space density calculation.")
   end

   # Map component names to columns
   component_names = data.header.real_component_names
   component_map = Dict(name => i for (i, name) in enumerate(component_names))

   x_variable = _resolve_alias(x_variable)
   y_variable = _resolve_alias(y_variable)

   if !haskey(component_map, x_variable) || !haskey(component_map, y_variable)
      error("Invalid variable name. Available: $(keys(component_map))")
   end

   x_index = component_map[x_variable]
   y_index = component_map[y_variable]

   x_data = rdata[:, x_index]
   y_data = rdata[:, y_index]

   # Using FHist.jl
   # FHist uses nbins keyword
   arg_bins = bins isa Int ? (bins, bins) : bins
   nx, ny = arg_bins
   h = Hist2D((x_data, y_data); nbins=arg_bins, overflow=false)
   # FHist might return extra bins, or OffsetArrays. We slice to expected range.
   # Assuming 1-based indexing for the main bins.
   H = collect(h.bincounts[1:nx, 1:ny])
   xedges = collect(h.binedges[1])
   yedges = collect(h.binedges[2])

   return H, xedges, yedges
end
