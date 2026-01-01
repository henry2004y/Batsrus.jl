# Plotting functionalities for AMReX particle data.

using FHist

export get_phase_space_density, plot_phase

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
   H = h.bincounts
   xedges = collect(h.binedges[1])
   yedges = collect(h.binedges[2])

   return H, xedges, yedges
end

"""
    plot_phase(data, x_var, y_var; bins=100, x_range=nothing, y_range=nothing, z_range=nothing, log_scale=true, ax=nothing, add_colorbar=true, kwargs...)

Plots the 2D phase space density for selected variables.
"""
function plot_phase(
      data::AMReXParticleData,
      x_variable::String,
      y_variable::String;
      bins::Union{Int, Tuple{Int, Int}} = 100,
      x_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      y_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      z_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
      log_scale::Bool = true,
      ax = nothing,
      add_colorbar::Bool = true,
      kwargs...
)
   # 1. Get phase space density
   H, xedges, yedges = get_phase_space_density(
      data, x_variable, y_variable; bins, x_range, y_range, z_range
   )

   # 2. Plot
   if isnothing(ax)
      ax = plt.gca()
   end

   # Handle log scale
   norm = nothing
   if log_scale
      if maximum(H) > 0
         vmin = minimum(H[H .> 0])
         vmax = maximum(H)
         norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
      end
   else
      norm = matplotlib.colors.Normalize()
   end

   # Prepare extent
   extent = [xedges[1], xedges[end], yedges[1], yedges[end]]

   im = ax.imshow(
      H',
      origin="lower",
      extent=extent,
      aspect="auto",
      norm=norm;
      kwargs...
   )

   if add_colorbar
      plt.colorbar(im, ax=ax, label="Count", pad=0.02)
   end

   xlabel(_get_axis_label(_resolve_alias(x_variable)))
   ylabel(_get_axis_label(_resolve_alias(y_variable)))

   return im
end
