# Plotting functionalities for AMReX particle data using Makie.

const _AXIS_LABEL_MAP = Dict(
   "velocity_x" => L"v_x",
   "velocity_y" => L"v_y",
   "velocity_z" => L"v_z"
)

_get_axis_label(variable_name::String) = get(_AXIS_LABEL_MAP, variable_name, variable_name)

import Batsrus: plot_phase!

"""
    plot_phase(data, x_variable, y_variable; bins=100, edges=nothing, x_range=nothing, y_range=nothing, z_range=nothing, log_scale=true, axis=(;), ax=nothing, figure=(;), add_colorbar=true, transform=nothing, plot_zero_lines=false, normalize=false, kwargs...)

Plots the 2D phase space density for selected variables using Makie.
This function wraps `get_phase_space_density` and delegates plotting to `Makie.plot`.

# Arguments

  - `data`: `AMReXParticle` data object.
  - `x_variable`: Name of the variable for the x-axis (e.g., "vx").
  - `y_variable`: Name of the variable for the y-axis (e.g., "vy").
  - `bins`: Number of bins for the histogram (default: 100).
  - `edges`: **Histogram binning edges**. Explicitly defines the bin edges for `x_variable` and `y_variable`. Overrides `bins`.
  - `x_range`, `y_range`, `z_range`: **Spatial selection ranges**. Only particles within these ranges in configuration space are included.
  - `log_scale`: Whether to use a logarithmic color scale (default: `true`).
  - `axis`: NamedTuple of keyword arguments passed to `Makie.Axis` (used when creating a new axis).
  - `ax`: Existing `Makie.Axis` to plot into.
  - `figure`: NamedTuple of keyword arguments passed to `Makie.Figure` (used when creating a new figure).
  - `add_colorbar`: Whether to add a colorbar to the plot (default: `true`). Only applies when creating a new figure.
  - `transform`: Optional function to transform the data before binning.
  - `plot_zero_lines`: Whether to draw dashed lines at x=0 and y=0 (default: `false`).
  - `normalize`: Whether to normalize the histogram to a probability density (default: `false`).
  - `kwargs`: Additional keyword arguments passed to `Makie.plot` (e.g., `colormap`).

TODO: support 1D/2D/3D phase space plotting.
"""
function Batsrus.plot_phase!(
      ax::Union{Makie.Axis, Nothing},
      data::AMReXParticle,
      x_variable::String,
      y_variable::String;
      bins::Union{Int, Tuple{Int, Int}} = 100,
      edges = nothing,
      x_range = nothing,
      y_range = nothing,
      z_range = nothing,
      log_scale::Bool = true,
      axis = (;),
      figure = (;),
      add_colorbar::Bool = true,
      transform::Union{Function, Nothing} = nothing,
      plot_zero_lines::Bool = false,
      normalize::Bool = false,
      kwargs...
)
   # Get phase space density (Hist2D object)
   h = Batsrus.get_phase_space_density(
      data, x_variable, y_variable; bins, edges, x_range, y_range, z_range, transform, normalize
   )

   plot_kwargs = Dict{Symbol, Any}(kwargs...)

   if log_scale
      #TODO log10 has some issues with the mock data
      plot_kwargs[:colorscale] = Makie.pseudolog10
   end

   xlabel = _get_axis_label(_resolve_alias(x_variable))
   ylabel = _get_axis_label(_resolve_alias(y_variable))

   if !isnothing(ax)
      # Plot into existing axis
      if ax.xlabel[] == "x" || ax.xlabel[] == ""
         ax.xlabel = xlabel
      end
      if ax.ylabel[] == "y" || ax.ylabel[] == ""
         ax.ylabel = ylabel
      end

      pl = Makie.heatmap!(ax, h; plot_kwargs...)

      if plot_zero_lines
         Makie.hlines!(ax, [0]; color = :gray, linestyle = :dash, linewidth = 1)
         Makie.vlines!(ax, [0]; color = :gray, linestyle = :dash, linewidth = 1)
      end

      return pl
   else
      # Create new figure and axis
      axis_kwargs = merge((xlabel = xlabel, ylabel = ylabel), axis)

      p = Makie.heatmap(h; axis = axis_kwargs, figure = figure, plot_kwargs...)

      if plot_zero_lines
         Makie.hlines!(p.axis, [0]; color = :gray, linestyle = :dash, linewidth = 1)
         Makie.vlines!(p.axis, [0]; color = :gray, linestyle = :dash, linewidth = 1)
      end

      if add_colorbar
         Makie.Colorbar(
            p.figure[1, 2], p.plot; label = normalize ? "Probability Density" : "Counts")
      end

      return p
   end
end
