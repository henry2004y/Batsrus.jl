# Type conversion from Batsrus to Makie

"""
Conversion for 1D plots
"""
function Makie.convert_arguments(
        P::Makie.PointBased, bd::BatsrusIDL, var::String;
        use_units = false
    )
    da = bd[var]

    if use_units && hasunit(bd)
        unitx = getunit(bd, bd.head.coord[1])
        unitw = getunit(bd, var)
        if unitx isa UnitfulBatsrus.Unitlike
            da = DimArray(parent(da), (dims(da, 1) * unitx,))
        end
        if unitw isa UnitfulBatsrus.Unitlike
            da = da .* unitw
        end
    end

    da_stripped = ustrip.(da)
    return (Point2f.(lookup(da_stripped, 1), parent(da_stripped)),)
end

"""
Conversion for 2D plots.
"""
function Makie.convert_arguments(
        P::Makie.GridBased, bd::BatsrusIDL, var::String;
        plotrange = [-Inf, Inf, -Inf, Inf], plotinterval = 0.1,
        use_units = false
    )
    x, y, w = interp2d(bd, var, plotrange, plotinterval)

    if use_units && hasunit(bd)
        unitx = getunit(bd, bd.head.coord[1])
        unity = getunit(bd, bd.head.coord[2])
        unitw = getunit(bd, var)

        if unitx isa UnitfulBatsrus.Unitlike
            x *= unitx
        end
        if unity isa UnitfulBatsrus.Unitlike
            y *= unity
        end
        if unitw isa UnitfulBatsrus.Unitlike
            w *= unitw
        end
    end

    # We still return a Tuple of arrays as expected by GridBased traits
    return (ustrip.(x), ustrip.(y), ustrip.(w)')
end
