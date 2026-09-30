# Sequence of figure plotting

export animate

"""
    animate(files::Vector{String}; kwargs...)

Save figures of colored contour from SWMF output files.
Uses DimensionalData selectors or interpolation for data extraction.

Unlike a single-snapshot plot, the color limits are determined by scanning **all**
snapshots when `vmin` or `vmax` is `Inf`, so the color scale stays consistent across
the whole animation.

# Keywords
- `var`: variable to plot with pcolormesh.
- `vmin::Real`: minimum plotting value (`Inf` = scan all files).
- `vmax::Real`: maximum plotting value (`Inf` = scan all files).
- `colorscale::Symbol`: color scale for the plot (`:linear`, `:log`, `:symlog`, `:twoslope`).
- `vcenter::Real`: center value of the color scale for `:twoslope`.
- `plotrange`: 2D plotting spatial range `[xmin, xmax, ymin, ymax]`.
- `plotinterval`: spatial sampling interval. `nothing` (default) auto-detects the
  finest grid resolution via [`finest_resolution`](@ref) — recommended for AMR and
  unstructured output. `Inf` keeps the original sampling.
- `innermask::Bool`: if true, mask the inner boundary (useful for generalized coordinates).
- `rbody`: inner body radius for the mask.
- `xlabel`, `ylabel`: strings for the x and y axis labels.
- `title`: title of the figure. Can be a `String`, or a function `(bd) -> String` that takes the `BatsrusIDL` object to generate time-dependent titles.
- `xlabel_kwargs`, `ylabel_kwargs`, `title_kwargs`, `clabel_kwargs`: keyword arguments for labels and title properties.
- `plot_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `pcolormesh`. If `nothing` (default), the colormap is selected automatically: `RdBu_r` with symmetric limits for bipolar data, `turbo` otherwise; override with `colormap`.
- `symmetric::Bool`: force symmetric color limits about zero (default: automatic for bipolar data).
- `colormap`: colormap overriding the automatic selection.
- `aspect`: aspect ratio passed to `set_aspect`; `nothing` keeps the default `equal`.
- `cbar_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `colorbar`. You can also specify a `label` here.
- `streamvars`: string with two variables separated by `;` for streamlines (e.g., `"bx;by"`).
- `stream_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `streamplot`.
- `outdir::String`: output directory for the image files.
- `overwrite::Bool`: if true, overwrite the existing image files.
- `fig_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `figure`.
- `savefig_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `savefig`.
- `cache_frames::Bool`: if true (default), frames extracted during the color-range scan are reused for rendering, avoiding a second read of every file.
- `show_progress::Bool`: if true, show progress bar.
"""
function Batsrus.animate(
        files::Vector{String}; var = "rho", vmin = -Inf, vmax = Inf, vcenter = 0.0,
        colorscale = :linear, plotrange = nothing, plotinterval = nothing,
        innermask = false, rbody = nothing,
        xlabel = nothing, ylabel = nothing, title = nothing,
        xlabel_kwargs = (;), ylabel_kwargs = (;), title_kwargs = (;),
        clabel_kwargs = (;),
        cbar_kwargs = (; orientation = "vertical", extend = "neither", pad = 0.005),
        streamvars = nothing,
        stream_kwargs = (; color = "white", density = 1),
        outdir = "figs/", overwrite = false,
        fig_kwargs = (; figsize = (8, 6), constrained_layout = true),
        savefig_kwargs = (; bbox_inches = "tight", dpi = 200),
        use_units = false,
        colormap = nothing, symmetric = nothing, aspect = nothing,
        plot_kwargs = nothing,
        cache_frames::Bool = true,
        show_progress = true,
        kwargs...
    )

    if !isdir(outdir)
        mkpath(outdir)
    end

    nfile = length(files)
    if nfile == 0
        @warn "No files found to process."
        return
    end

    var = var isa AbstractString ? String(var) : var

    # Resolve plotinterval: `nothing` auto-detects the finest cell size from the
    # first snapshot (SMR assumption — the mesh is static, so one file is
    # representative of the entire run).
    dx_fine = if isnothing(plotinterval)
        finest = finest_resolution(load(files[1]))
        if show_progress
            @info "Auto-detected finest grid resolution: plotinterval = $(round(finest, sigdigits = 5))"
        end
        Float64(finest)
    elseif isinf(plotinterval)
        Inf32  # legacy behavior: no resampling
    else
        Float64(plotinterval)
    end

    # Extract a plottable 2D frame (x coords, y coords, data) from `bd`.
    # Always uses the current value of dx_fine (updated by auto-detection above).
    function _extract(bd, v)
        if bd.head.gencoord
            extents = isnothing(plotrange) ? [-Inf, Inf, -Inf, Inf] : plotrange
            x_coords, y_coords, data = interp2d(bd, v, extents, dx_fine; innermask, rbody)
        else
            data = bd[v]
            if !isnothing(plotrange)
                sel1, sel2 = plotrange[1] .. plotrange[2], plotrange[3] .. plotrange[4]
                data = data[sel1, sel2]
            end
            x_coords = dims(data, 1).val
            y_coords = dims(data, 2).val
            # Ensure uniform spacing for Matplotlib with Float64 precision
            x_coords = range(Float64(x_coords[1]), Float64(x_coords[end]), length(x_coords))
            y_coords = range(Float64(y_coords[1]), Float64(y_coords[end]), length(y_coords))
            data = collect(data')
        end

        if use_units && hasunit(bd)
            unitx = getunit(bd, bd.head.coord[1])
            unity = getunit(bd, bd.head.coord[2])
            unitw = getunit(bd, v)
            if unitx isa UnitfulBatsrus.Unitlike
                x_coords = x_coords .* unitx
            end
            if unity isa UnitfulBatsrus.Unitlike
                y_coords = y_coords .* unity
            end
            if unitw isa UnitfulBatsrus.Unitlike
                data = data .* unitw
            end
        end

        return x_coords, y_coords, data
    end

    # Scan all snapshots to determine a consistent global color range.
    frames = Vector{Any}(undef, nfile)
    if isinf(vmin) || isinf(vmax)
        if show_progress
            @info "Scanning $(nfile) snapshots to determine range for '$var'..."
        end
        global_min = Inf
        global_max = -Inf
        p = Progress(nfile; dt = 2, enabled = show_progress)
        for (i, file) in enumerate(files)
            bd = load(file)
            frame = _extract(bd, var)
            cache_frames && (frames[i] = frame)
            cur_data = frame[end]
            valid_data = if colorscale == :log
                filter(x -> isfinite(x) && x > 0, cur_data)
            else
                filter(isfinite, cur_data)
            end
            if !isempty(valid_data)
                global_min = min(global_min, minimum(valid_data))
                global_max = max(global_max, maximum(valid_data))
            end
            next!(p)
        end

        if isinf(vmin)
            vmin = isinf(global_min) ? (colorscale == :log ? 1.0e-6 : 0.0) : global_min
        end
        if isinf(vmax)
            vmax = isinf(global_max) ? 1.0 : global_max
        end
    end

    if vmin == vmax
        vmin -= 0.1 * abs(vmin) + 1.0e-4
        vmax += 0.1 * abs(vmax) + 1.0e-4
    end

    # Auto-select colormap and optionally symmetrize limits (if not overridden by caller)
    if isnothing(plot_kwargs)
        is_bipolar = vmin < 0 && vmax > 0
        do_symmetric = isnothing(symmetric) ? is_bipolar : symmetric

        selected_cmap = if !isnothing(colormap)
            colormap
        elseif do_symmetric && colorscale !== :log
            PyPlot.matplotlib.cm.RdBu_r
        else
            PyPlot.matplotlib.cm.turbo
        end

        if do_symmetric && colorscale !== :log
            vlim = max(abs(vmin), abs(vmax))
            vmin, vmax = -vlim, vlim
        end

        plot_kwargs = (; cmap = selected_cmap)
    end

    norm = set_colorbar(colorscale, vmin, vmax, [1.0]; vcenter)

    fig = figure(; fig_kwargs...)
    ax = plt.axes()

    c = nothing
    cb = nothing
    st = nothing
    p = Progress(nfile; dt = 4, enabled = show_progress)

    for (i, file) in enumerate(files)
        # Generate output name
        base, ext = splitext(basename(file))
        outname = joinpath(outdir, base * ".png")

        if !overwrite && isfile(outname)
            next!(p)
            continue
        end

        bd = load(file)
        x_coords, y_coords, data = if !isnothing(frames[i])
            frames[i]
        else
            _extract(bd, var)
        end

        if isnothing(c)
            # Initialization
            # Plotting
            c = ax.pcolormesh(x_coords, y_coords, data; norm, plot_kwargs...)
            # Labels and aspect ratio
            x_label = isnothing(xlabel) ? L"X [$R_\mathrm{E}$]" : xlabel
            ax.set_xlabel(x_label; xlabel_kwargs...)

            if isnothing(ylabel)
                fname = basename(files[1])
                y_label = startswith(fname, "y") ? L"Z [$R_\mathrm{E}$]" : L"Y [$R_\mathrm{E}$]"
            else
                y_label = ylabel
            end
            ax.set_ylabel(y_label; ylabel_kwargs...)

            if isnothing(aspect)
                ax.set_aspect("equal", adjustable = "box", anchor = "C")
            else
                ax.set_aspect(aspect)
            end
            ax.set_xlim(x_coords[1], x_coords[end])
            ax.set_ylim(y_coords[1], y_coords[end])

            if isnothing(cb)
                cb = colorbar(c; ax, cbar_kwargs...)
                c_label = hasproperty(cbar_kwargs, :label) ? cbar_kwargs.label : var
                cb.ax.set_ylabel(c_label; clabel_kwargs...)
            end
        else
            # Optimization: only update the array data
            c.set_array(data)
        end

        # Update Title
        title_str = if isnothing(title)
            @sprintf "t = %.1f s" bd.head.time
        elseif title isa Function
            title(bd)
        else
            title
        end
        ax.set_title(title_str; title_kwargs...)

        if !isnothing(streamvars)
            vars = split(streamvars, ";")
            if length(vars) == 2
                try
                    if bd.head.gencoord
                        extents = isnothing(plotrange) ? [-Inf, Inf, -Inf, Inf] : plotrange
                        xi, yi, (v1, v2) = interp2d(
                            bd, [String(vars[1]), String(vars[2])], extents, dx_fine;
                            innermask, rbody
                        )
                    else
                        v1 = bd[vars[1]]
                        v2 = bd[vars[2]]
                        if !isnothing(plotrange)
                            sel1, sel2 = plotrange[1] .. plotrange[2], plotrange[3] .. plotrange[4]
                            v1 = v1[sel1, sel2]
                            v2 = v2[sel1, sel2]
                        end
                        v1, v2 = collect(v1'), collect(v2')
                        xi, yi = x_coords, y_coords
                    end
                    # Assuming same grid for streamlines
                    st = ax.streamplot(xi, yi, v1, v2; stream_kwargs...)
                catch err
                    @debug "Streamplot skipped for $(vars[1]);$(vars[2]): $(err)"
                end
            end
        end

        savefig(outname; savefig_kwargs...)

        if !isnothing(st) # Clean up streamlines for next iteration
            st.lines.remove()
            for art in ax.get_children()
                if PyPlot.PyCall.pybuiltin(:isinstance)(
                        art,
                        PyPlot.matplotlib.patches.FancyArrowPatch
                    )
                    art.remove()
                end
            end
            st = nothing
        end
        next!(p)
    end

    empty!(frames) # release cached frames
    PyPlot.plt.close(fig)
    return
end
