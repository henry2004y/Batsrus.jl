# Sequence of figure plotting

export animate

"""
    animate(files::Vector{String}; kwargs...)

Save figures of colored contour from SWMF output files.
Uses DimensionalData selectors or interpolation for data extraction.

# Keywords
- `var::String`: variable to plot with pcolormesh.
- `vmin::Real`: minimum plotting value.
- `vmax::Real`: maximum plotting value.
- `plotrange`: 2D plotting spatial range `[xmin, xmax, ymin, ymax]`.
- `plotinterval`: spatial sampling interval.
- `innermask::Bool`: if true, mask the inner boundary (useful for generalized coordinates).
- `rbody`: inner body radius for the mask.
- `xlabel`, `ylabel`: strings for the x and y axis labels.
- `title`: title of the figure. Can be a `String`, or a function `(bd) -> String` that takes the `BatsrusIDL` object to generate time-dependent titles.
- `plot_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `pcolormesh`.
- `cbar_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `colorbar`. You can also specify a `label` here.
- `streamvars`: string with two variables separated by `;` for streamlines (e.g., `"bx;by"`).
- `stream_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `streamplot`.
- `outdir::String`: output directory for the image files.
- `overwrite::Bool`: if true, overwrite the existing image files.
- `fig_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `figure`.
- `savefig_kwargs`: NamedTuple or Dict of keyword arguments passed to Matplotlib's `savefig`.
"""
function Batsrus.animate(
        files::Vector{String}; var = "rho", vmin = -Inf, vmax = Inf,
        plotrange = nothing, plotinterval = Inf,
        innermask = false, rbody = nothing,
        xlabel = nothing, ylabel = nothing, title = nothing,
        plot_kwargs = (; cmap = PyPlot.matplotlib.cm.turbo),
        cbar_kwargs = (; orientation = "vertical", extend = "neither", pad = 0.005),
        streamvars = nothing,
        stream_kwargs = (; color = "white", density = 1),
        outdir = "figs/", overwrite = false,
        fig_kwargs = (; figsize = (8, 6), constrained_layout = true),
        savefig_kwargs = (; bbox_inches = "tight", dpi = 200),
        use_units = false
    )

    if !isdir(outdir)
        mkpath(outdir)
    end

    nfile = length(files)
    if nfile == 0
        @warn "No files found to process."
        return
    end

    # First file to setup parameters if vmin/vmax are Inf
    bd = load(files[1])
    if bd.head.gencoord
        range = isnothing(plotrange) ? [-Inf, Inf, -Inf, Inf] : plotrange
        _, _, data = interp2d(bd, var, range, plotinterval; innermask, rbody)
    else
        data = bd[var]
        if !isnothing(plotrange)
            data = data[
                dims(data, 1) => plotrange[1] .. plotrange[2],
                dims(data, 2) => plotrange[3] .. plotrange[4],
            ]
        end
    end

    if isinf(vmin) || isinf(vmax)
        vmin = isinf(vmin) ? minimum(data) : vmin
        vmax = isinf(vmax) ? maximum(data) : vmax
    end

    norm = PyPlot.matplotlib.colors.Normalize(vmin, vmax)

    fig = figure(; fig_kwargs...)
    ax = plt.axes()

    c = nothing
    cb = nothing
    st = nothing
    for (i, file) in enumerate(files)
        @info "Processing $i out of $nfile: $file"

        # Generate output name
        base, ext = splitext(basename(file))
        outname = joinpath(outdir, base * ".png")

        if !overwrite && isfile(outname)
            @info "Skipping existing file: $outname"
            continue
        end

        bd = load(file)
        if bd.head.gencoord
            range = isnothing(plotrange) ? [-Inf, Inf, -Inf, Inf] : plotrange
            x_coords, y_coords, data =
                interp2d(bd, var, range, plotinterval; innermask, rbody)
        else
            data = bd[var]
            # Apply plotrange if provided (using DimensionalData selectors)
            if !isnothing(plotrange)
                data = data[
                    dims(data, 1) => plotrange[1] .. plotrange[2],
                    dims(data, 2) => plotrange[3] .. plotrange[4],
                ]
            end
            # Get coordinates
            x_coords = dims(data, 1).val
            y_coords = dims(data, 2).val
            data = data'
        end

        if use_units && hasunit(bd)
            unitx = getunit(bd, bd.head.coord[1])
            unity = getunit(bd, bd.head.coord[2])
            unitw = getunit(bd, var)
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

        if isnothing(c)
            # Initialization
            # Plotting
            c = ax.pcolormesh(x_coords, y_coords, data; norm, plot_kwargs...)
            # Labels and aspect ratio
            x_label = isnothing(xlabel) ? L"X [$R_\mathrm{E}$]" : xlabel
            ax.set_xlabel(x_label)

            if isnothing(ylabel)
                fname = basename(files[1])
                y_label = startswith(fname, "y") ? L"Z [$R_\mathrm{E}$]" : L"Y [$R_\mathrm{E}$]"
            else
                y_label = ylabel
            end
            ax.set_ylabel(y_label)

            ax.set_aspect("equal", adjustable = "box", anchor = "C")
            ax.set_xlim(x_coords[1], x_coords[end])
            ax.set_ylim(y_coords[1], y_coords[end])

            if isnothing(cb)
                cb = colorbar(c; ax, cbar_kwargs...)
                c_label = hasproperty(cbar_kwargs, :label) ? cbar_kwargs.label : var
                cb.ax.set_ylabel(c_label)
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
        ax.set_title(title_str)

        if !isnothing(streamvars)
            vars = split(streamvars, ";")
            if length(vars) == 2
                if bd.head.gencoord
                    range = isnothing(plotrange) ? [-Inf, Inf, -Inf, Inf] : plotrange
                    xi, yi, v1 =
                        interp2d(bd, String(vars[1]), range, plotinterval; innermask, rbody)
                    _, _, v2 =
                        interp2d(bd, String(vars[2]), range, plotinterval; innermask, rbody)
                else
                    v1 = bd[vars[1]]
                    v2 = bd[vars[2]]
                    if !isnothing(plotrange)
                        v1 = v1[
                            dims(v1, 1) => plotrange[1] .. plotrange[2],
                            dims(v1, 2) => plotrange[3] .. plotrange[4],
                        ]
                        v2 = v2[
                            dims(v2, 1) => plotrange[1] .. plotrange[2],
                            dims(v2, 2) => plotrange[3] .. plotrange[4],
                        ]
                    end
                    v1, v2 = collect(v1'), collect(v2')
                    xi, yi = x_coords, y_coords
                end
                # Assuming same grid for streamlines
                st = ax.streamplot(xi, yi, v1, v2; stream_kwargs...)
            end
        end

        savefig(outname; savefig_kwargs...)
        @info "Saved $outname"

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
    end

    return PyPlot.plt.close(fig)
end
