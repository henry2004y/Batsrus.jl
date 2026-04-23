# Sequence of figure plotting with Makie

export animate

"""
    animate(files::Vector{String}; kwargs...)

Save figures of colored contour from SWMF output files using Makie.
Uses DimensionalData selectors or interpolation for data extraction.

# Keywords
- `var::String`: variable to plot with heatmap.
- `vmin::Real`: minimum plotting value.
- `vmax::Real`: maximum plotting value.
- `plotrange`: 2D plotting spatial range `[xmin, xmax, ymin, ymax]`.
- `plotinterval`: spatial sampling interval.
- `innermask::Bool`: if true, mask the inner boundary (useful for generalized coordinates).
- `rbody`: inner body radius for the mask.
- `xlabel`, `ylabel`: strings for the x and y axis labels.
- `title`: title of the figure. Can be a `String`, or a function `(bd) -> String` that
takes the `BatsrusIDL` object to generate time-dependent titles.
- `colormap`: colormap passed to Makie's `heatmap`.
- `streamvars`: string with two variables separated by `;` for streamlines (e.g., `"bx;by"`).
- `stream_kwargs`: NamedTuple or Dict of keyword arguments passed to `streamlines!`.
- `outdir::String`: output directory for the image files.
- `overwrite::Bool`: if true, overwrite the existing image files.
- `fig_kwargs`: NamedTuple or Dict of keyword arguments passed to Makie's `Figure`.
- `savefig_kwargs`: NamedTuple or Dict of keyword arguments passed to Makie's `save`.
- `colgap`: gap between the axis and the colorbar.
"""
function Batsrus.animate(
        files::Vector{String}; var = "rho", vmin = -Inf, vmax = Inf,
        plotrange = nothing, plotinterval = 0.5,
        innermask = false, rbody = nothing,
        xlabel = nothing, ylabel = nothing, title = nothing,
        colormap = :turbo,
        streamvars = nothing,
        stream_kwargs = (color = :white, linewidth = 1.0, with_arrows = true),
        outdir = "figs/", overwrite = false,
        fig_kwargs = (size = (800, 600),),
        savefig_kwargs = Dict(),
        use_units = false,
        colgap = 10
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
        _, _, data = interp2d(
            bd, var, range, plotinterval; innermask, rbody,
            useMatplotlib = false
        )
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

    fig = Figure(; fig_kwargs...)

    # Determine labels
    x_label = isnothing(xlabel) ? L"X [$R_\mathrm{E}$]" : xlabel
    if isnothing(ylabel)
        fname = basename(files[1])
        y_label = startswith(fname, "y") ? L"Z [$R_\mathrm{E}$]" : L"Y [$R_\mathrm{E}$]"
    else
        y_label = ylabel
    end

    ax = Axis(fig[1, 1], xlabel = x_label, ylabel = y_label, aspect = DataAspect())

    # Initialize plot objects
    hm = nothing
    cb = nothing
    st = nothing
    inner_mask_poly = nothing

    for (i, file) in enumerate(files)
        @info "Processing $i out of $nfile: $file"

        # Generate output name
        base, _ = splitext(basename(file))
        outname = joinpath(outdir, base * ".png")

        if !overwrite && isfile(outname)
            @info "Skipping existing file: $outname"
            continue
        end

        bd = load(file)
        if bd.head.gencoord
            range = isnothing(plotrange) ? [-Inf, Inf, -Inf, Inf] : plotrange
            x_coords, y_coords, data =
                interp2d(
                bd, var, range, plotinterval; innermask, rbody,
                useMatplotlib = false
            )
        else
            data = bd[var]
            if !isnothing(plotrange)
                data = data[
                    dims(data, 1) => plotrange[1] .. plotrange[2],
                    dims(data, 2) => plotrange[3] .. plotrange[4],
                ]
            end
            x_coords = dims(data, 1).val
            y_coords = dims(data, 2).val
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

        if isnothing(hm)
            hm = heatmap!(
                ax, x_coords, y_coords, data';
                colormap, colorrange = (vmin, vmax)
            )
            cb = Colorbar(fig[1, 2], hm, label = var)
            # Adjust gap between axis and colorbar
            colgap!(fig.layout, colgap)
        else
            # For simplicity in this implementation, we recreate the plot elements
            # that might change significantly, like streamlines.
            # Heatmap data update:
            hm[3] = data'
        end

        # Update Title
        ax.title = isnothing(title) ? @sprintf("t = %.1f s", bd.head.time) :
            (title isa Function ? title(bd) : title)

        # Streamlines
        if !isnothing(streamvars)
            if !isnothing(st)
                delete!(ax, st)
            end
            vars = split(streamvars, ";")
            if length(vars) == 2
                if bd.head.gencoord
                    _, _, v1 = interp2d(
                        bd, String(vars[1]), range, plotinterval;
                        innermask, rbody, useMatplotlib = false
                    )
                    _, _, v2 = interp2d(
                        bd, String(vars[2]), range, plotinterval;
                        innermask, rbody, useMatplotlib = false
                    )
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
                    v1, v2 = v1.val, v2.val
                end
                str = evenstream(x_coords, y_coords, v1', v2')
                st = streamlines!(ax, str; stream_kwargs...)
            end
        end

        # Inner Boundary Mask
        if innermask && !isnothing(rbody)
            if isnothing(inner_mask_poly)
                inner_mask_poly = poly!(ax, Circle(Point2f(0, 0), rbody), color = :white)
            end
        end

        save(outname, fig; savefig_kwargs...)
        @info "Saved $outname"
    end

    return
end
