# Plotting functionalities for AMReX particle data.

export plot_phase

const _AXIS_LABEL_MAP = Dict(
   "velocity_x" => "\$v_x\$",
   "velocity_y" => "\$v_y\$",
   "velocity_z" => "\$v_z\$"
)

_get_axis_label(variable_name::String) = get(_AXIS_LABEL_MAP, variable_name, variable_name)

import Batsrus: plot_phase!

"""
    plot_phase(data, x_variable, y_variable; bins=100, edges=nothing, x_range=nothing, y_range=nothing, z_range=nothing, log_scale=true, ax=nothing, add_colorbar=true, transform=nothing, plot_zero_lines=false, normalize=false, kwargs...)

Plots the 2D phase space density for selected variables.

# Arguments

  - `data`: `AMReXParticle` data object.
  - `x_variable`: Name of the variable for the x-axis (e.g., "vx").
  - `y_variable`: Name of the variable for the y-axis (e.g., "vy").
  - `bins`: Number of bins for the histogram (default: 100).
  - `edges`: **Histogram binning edges**. Explicitly defines the bin edges for `x_variable` and `y_variable`. Overrides `bins`.
  - `x_range`, `y_range`, `z_range`: **Spatial selection ranges**. Only particles within these ranges in configuration space are included.
  - `log_scale`: Whether to use a logarithmic color scale (default: `true`).
  - `ax`: PyPlot axis to plot on. If `nothing`, uses the current axis.
  - `add_colorbar`: Whether to add a colorbar to the plot (default: `true`).
  - `transform`: Optional function to transform the data before binning.
  - `plot_zero_lines`: Whether to draw dashed lines at x=0 and y=0 (default: `false`).
  - `normalize`: Whether to normalize the histogram to a probability density (default: `false`).
  - `vmin`, `vmax`: **Histogram range**. Explicitly defines the range for the histogram.
  - `kwargs`: Additional keyword arguments passed to `imshow` (e.g., `cmap`).
"""
function Batsrus.plot_phase!(
      ax::Union{PyPlot.PyObject, Nothing},
      data::AMReXParticle,
      x_variable::String,
      y_variable::String;
      bins::Union{Int, Tuple{Int, Int}} = 100,
      edges = nothing,
      x_range = nothing,
      y_range = nothing,
      z_range = nothing,
      log_scale::Bool = true,
      add_colorbar::Bool = true,
      transform::Union{Function, Nothing} = nothing,
      plot_zero_lines::Bool = false,
      normalize::Bool = false,
      vmin = nothing,
      vmax = nothing,
      kwargs...
)
   # Get phase space density
   h = get_phase_space_density(
      data, x_variable, y_variable; bins, edges, x_range, y_range, z_range, transform, normalize
   )

   H = h.bincounts
   xedges = h.binedges[1]
   yedges = h.binedges[2]

   # Plot
   if isnothing(ax)
      ax = plt.gca()
   end

   # Handle log scale
   norm = nothing
   if log_scale
      if maximum(H) > 0
         if isnothing(vmin)
            vmin = minimum(H[H .> 0])
         end
         if isnothing(vmax)
            vmax = maximum(H)
         end
         norm = PyPlot.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax)
      end
   else
      norm = PyPlot.matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
   end

   extent = [xedges[1], xedges[end], yedges[1], yedges[end]]

   im = ax.imshow(
      H';
      origin = "lower",
      extent,
      aspect = "auto",
      norm,
      kwargs...
   )

   if plot_zero_lines
      ax.axhline(0, color = "gray", linestyle = "--", linewidth = 1, alpha = 0.5)
      ax.axvline(0, color = "gray", linestyle = "--", linewidth = 1, alpha = 0.5)
   end

   if add_colorbar
      cbar_label = normalize ? "Probability Density" : "PSD"
      plt.colorbar(im; ax, label = cbar_label, pad = 0.02)
   end

   xlabel(_get_axis_label(_resolve_alias(x_variable)))
   ylabel(_get_axis_label(_resolve_alias(y_variable)))

   return im
end
