# Plotting functionalities for AMReX particle data.


export plot_phase


const _AXIS_LABEL_MAP = Dict(
   "velocity_x" => "\$v_x\$",
   "velocity_y" => "\$v_y\$",
   "velocity_z" => "\$v_z\$"
)



_get_axis_label(variable_name::String) = get(_AXIS_LABEL_MAP, variable_name, variable_name)



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
   # Get phase space density
   H, xedges, yedges = get_phase_space_density(
      data, x_variable, y_variable; bins, x_range, y_range, z_range
   )

   # Plot
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
