# Utility functions for plotting and analyzing.

"""
    interp2d(bd::BatsrusIDL, var, plotrange=[-Inf, Inf, -Inf, Inf],
     	plotinterval=Inf; kwargs...)

Return 2D interpolated slices of data `var` from `bd`. If `plotrange` is not set, output
data resolution is the same as the original.

`var` can be an `AbstractString` for file variables or a `Symbol` for derived quantities
(e.g. `:b`, `:anisotropy`). Using a `Symbol` is recommended for performance-critical
calls because the returned array type is fully resolved at compile time.

# Keyword Arguments

  - `innermask=false`: Whether to mask the inner boundary with NaN.
  - `rbody=nothing`: Radius of the inner mask. If not set, it will be read from the header.
  - `useMatplotlib=true`: Whether to Matplotlib (faster) or NaturalNeighbours for scattered
    interpolation. If true, a linear interpolation is performed on a constructed triangle mesh.
"""
function interp2d end

function interp2d(
        bd::BatsrusIDLUnstructured{2, TV, TX, TW},
        var::Union{AbstractString, Symbol},
        plotrangeIn::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32;
        innermask::Bool = false, rbody::Union{Nothing, Real} = nothing,
        useMatplotlib::Bool = true
    ) where {TV, TX, TW}
    W_raw = getvar(bd, var)
    xi, yi, Wis = _interp2d_unstructured(
        bd, [W_raw], plotrangeIn, plotinterval; innermask, rbody, useMatplotlib
    )
    return xi, yi, Wis[1]
end

function interp2d(
        bd::BatsrusIDLUnstructured{2, TV, TX, TW},
        vars::AbstractVector{<:Union{AbstractString, Symbol}},
        plotrangeIn::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32;
        innermask::Bool = false, rbody::Union{Nothing, Real} = nothing,
        useMatplotlib::Bool = true
    ) where {TV, TX, TW}
    Ws_raw = [getvar(bd, v) for v in vars]
    return _interp2d_unstructured(
        bd, Ws_raw, plotrangeIn, plotinterval; innermask, rbody, useMatplotlib
    )
end

function _interp2d_unstructured(
        bd::BatsrusIDLUnstructured{2, TV, TX, TW}, Ws_raw::AbstractVector,
        plotrange::Vector, plotinterval::Real;
        innermask::Bool, rbody::Union{Nothing, Real}, useMatplotlib::Bool
    ) where {TV, TX, TW}
    x = bd.x
    X, Y = eachslice(x, dims = 3)
    X, Y = vec(X), vec(Y)
    Ws = [vec(W_raw) for W_raw in Ws_raw]

    adjust_plotrange!(plotrange, extrema(X), extrema(Y))
    if isinf(plotinterval)
        plotinterval = (plotrange[2] - plotrange[1]) / size(X, 1)
    end
    xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
    yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)

    Wis = [
        if useMatplotlib
                _triangulate_matplotlib(X, Y, W, xi, yi)
        else
                _, _, Wi_ = interpolate2d_generalized_coords(X, Y, W, plotrange, plotinterval)
                Wi_
        end for W in Ws
    ]

    for Wi in Wis
        _mask_inner_boundary!(Wi, xi, yi, bd, innermask, rbody)
    end

    return xi, yi, Wis
end

function interp2d(
        bd::BatsrusIDLStructured{2, TV, TX, TW},
        var::Union{AbstractString, Symbol},
        plotrangeIn::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32;
        innermask::Bool = false, rbody::Union{Nothing, Real} = nothing,
        useMatplotlib::Bool = true
    ) where {TV, TX, TW}
    W_raw = getvar(bd, var)
    return _interp2d_structured(
        bd, W_raw, plotrangeIn, plotinterval; innermask, rbody
    )
end

function interp2d(
        bd::BatsrusIDLStructured{2, TV, TX, TW},
        vars::AbstractVector{<:Union{AbstractString, Symbol}},
        plotrangeIn::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32;
        innermask::Bool = false, rbody::Union{Nothing, Real} = nothing
    ) where {TV, TX, TW}
    Ws_raw = [getvar(bd, v) for v in vars]
    return _interp2d_structured(
        bd, Ws_raw, plotrangeIn, plotinterval; innermask, rbody
    )
end

function _interp2d_structured(
        bd::BatsrusIDLStructured{2, TV, TX, TW}, Ws_raw::AbstractVector,
        plotrange::Vector, plotinterval::Real;
        innermask::Bool, rbody::Union{Nothing, Real}
    ) where {TV, TX, TW}
    xrange, yrange = get_range(bd)
    xi, yi, Wis = if all(isinf.(plotrange))
        xi_ = xrange
        yi_ = yrange
        Wis_ = [collect(W_raw)' for W_raw in Ws_raw]
        xi_, yi_, Wis_
    else
        adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

        xi_ = if isinf(plotinterval)
            range(plotrange[1], stop = plotrange[2], length = length(xrange))
        else
            range(plotrange[1], stop = plotrange[2], step = plotinterval)
        end
        yi_ = if isinf(plotinterval)
            range(plotrange[3], stop = plotrange[4], length = length(yrange))
        else
            range(plotrange[3], stop = plotrange[4], step = plotinterval)
        end
        Xf = repeat(TV.(xi_), inner = length(yi_))
        Yf = repeat(TV.(yi_), outer = length(xi_))
        Wis_ = [
            begin
                    itp = cubic_interp((xrange, yrange), parent(W))
                    Wif = Vector{TV}(undef, length(Xf))
                    itp(Wif, (Xf, Yf))
                    reshape(Wif, length(yi_), length(xi_))
                end for W in Ws_raw
        ]
        xi_, yi_, Wis_
    end

    for Wi in Wis
        _mask_inner_boundary!(Wi, xi, yi, bd, innermask, rbody)
    end

    return xi, yi, Wis
end

function _interp2d_structured(
        bd::BatsrusIDLStructured{2, TV, TX, TW}, W_raw,
        plotrange::Vector, plotinterval::Real;
        innermask::Bool, rbody::Union{Nothing, Real}
    ) where {TV, TX, TW}
    xrange, yrange = get_range(bd)
    xi, yi, Wi = if all(isinf.(plotrange))
        xi_ = xrange
        yi_ = yrange
        Wi_ = collect(W_raw)'
        xi_, yi_, Wi_
    else
        adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

        xi_ = if isinf(plotinterval)
            range(plotrange[1], stop = plotrange[2], length = length(xrange))
        else
            range(plotrange[1], stop = plotrange[2], step = plotinterval)
        end
        yi_ = if isinf(plotinterval)
            range(plotrange[3], stop = plotrange[4], length = length(yrange))
        else
            range(plotrange[3], stop = plotrange[4], step = plotinterval)
        end
        itp = cubic_interp((xrange, yrange), parent(W_raw))
        Xf = repeat(TV.(xi_), inner = length(yi_))
        Yf = repeat(TV.(yi_), outer = length(xi_))
        Wif = Vector{TV}(undef, length(Xf))
        itp(Wif, (Xf, Yf))
        Wi_ = reshape(Wif, length(yi_), length(xi_))
        xi_, yi_, Wi_
    end

    _mask_inner_boundary!(Wi, xi, yi, bd, innermask, rbody)

    return xi, yi, Wi
end

function _mask_inner_boundary!(
        Wi, xi, yi, bd::BatsrusIDL, innermask::Bool, rbody::Union{Nothing, Real}
    )
    !innermask && return
    if isnothing(rbody)
        varIndex_ = findlast(x -> x == "rbody", bd.head.param)
        if !isnothing(varIndex_)
            rbody = bd.head.eqpar[varIndex_]
        else
            @info "rbody not found in file header parameters; use default rbody=1.0"
            rbody = 1.0
        end
    end

    @inbounds @simd for i in CartesianIndices(Wi)
        # Wi[j, i] corresponds to yi[j] and xi[i]
        # i[1] is row index (j), i[2] is column index (i)
        if xi[i[2]]^2 + yi[i[1]]^2 < rbody^2
            Wi[i] = NaN
        end
    end

    return
end

"""
Return the axis range for 2D outputs. See [`interp2d`](@ref).
"""
function meshgrid(
        bd::BatsrusIDLUnstructured,
        plotrange::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32
    )
    x = bd.x

    X, Y = eachslice(x, dims = 3)
    X, Y = vec(X), vec(Y)

    adjust_plotrange!(plotrange, extrema(X), extrema(Y))
    # Set a heuristic value if not set
    if isinf(plotinterval)
        plotinterval = (plotrange[2] - plotrange[1]) / size(X, 1)
    end
    xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
    yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)

    return xi, yi
end

function meshgrid(
        bd::BatsrusIDLStructured,
        plotrange::Vector = [-Inf32, Inf32, -Inf32, Inf32], plotinterval::Real = Inf32
    )
    x = bd.x

    xrange, yrange = get_range(bd)
    if all(isinf.(plotrange))
        xi, yi = xrange, yrange
    else
        adjust_plotrange!(plotrange, (xrange[1], xrange[end]), (yrange[1], yrange[end]))

        if isinf(plotinterval)
            xi = range(plotrange[1], stop = plotrange[2], step = step(xrange))
            yi = range(plotrange[3], stop = plotrange[4], step = step(yrange))
        else
            xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
            yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
        end
    end

    return xi, yi
end

"""
Find variable index in the BATSRUS data.
"""
function findindex(bd::BatsrusIDL, var::AbstractString)
    for (i, name) in enumerate(bd.head.wname)
        if length(name) == length(var)
            match = true
            for (c1, c2) in zip(name, var)
                if lowercase(c1) != lowercase(c2)
                    match = false
                    break
                end
            end
            match && return i
        end
    end
    error("$(var) not found in file header variables!")
end

"""
Generating consistent 2D arrays for passing to plotting functions.
"""
function meshgrid(x, y)
    X = [x for _ in y, x in x]
    Y = [y for y in y, _ in x]

    return X, Y
end

@inline hasunit(bd::BatsrusIDL) = startswith(bd.head.headline, "normalized") ? false : true

"""
Adjust 2D plot ranges.
"""
function adjust_plotrange!(plotrange, xlimit, ylimit)
    plotrange[1] = ifelse(isinf(plotrange[1]), xlimit[1], plotrange[1])
    plotrange[2] = ifelse(isinf(plotrange[2]), xlimit[2], plotrange[2])
    plotrange[3] = ifelse(isinf(plotrange[3]), ylimit[1], plotrange[3])
    plotrange[4] = ifelse(isinf(plotrange[4]), ylimit[2], plotrange[4])

    return
end

"""
Perform Triangle interpolation of 2D data `W` on grid `X`, `Y`.
"""
function interpolate2d_generalized_coords(
        X::T, Y::T, W::T,
        plotrange::Vector{<:AbstractFloat}, plotinterval::Real
    ) where {T <: AbstractVector}
    xi = range(plotrange[1], stop = plotrange[2], step = plotinterval)
    yi = range(plotrange[3], stop = plotrange[4], step = plotinterval)
    itp = NN.interpolate(X, Y, W)
    #TODO Interpolate as a whole using predicates with multithreading
    Wi = [itp(x, y; method = NN.Triangle()) for y in yi, x in xi]::Matrix{eltype(W)}

    return xi, yi, Wi
end

"""
    interp1d(bd::BatsrusIDLStructured, var, loc::AbstractVector{<:AbstractFloat})

Interpolate `var` at spatial point `loc` in `bd`. `var` can be an `AbstractString` or
a `Symbol` for derived quantities.
"""
function interp1d(
        bd::BatsrusIDLStructured{2, TV, TX, TW},
        var::Union{AbstractString, Symbol},
        loc::AbstractVector{<:AbstractFloat}
    ) where {TV, TX, TW}
    v = getvar(bd, var)
    return _interp1d_point(bd, v, loc)
end

function _interp1d_point(
        bd::BatsrusIDLStructured{2, TV, TX, TW},
        v, loc::AbstractVector{<:AbstractFloat}
    ) where {TV, TX, TW}
    xrange, yrange = get_range(bd)
    itp = linear_interp((xrange, yrange), parent(v))

    return itp(Tuple(loc))
end

"""
    interp1d(bd::BatsrusIDLStructured, var, point1::Vector, point2::Vector)

Interpolate `var` along a line from `point1` to `point2` in `bd`. `var` can be an
`AbstractString` or a `Symbol` for derived quantities.
"""
function interp1d(
        bd::BatsrusIDLStructured{2, TV, TX, TW},
        var::Union{AbstractString, Symbol}, point1, point2
    ) where {TV, TX, TW}
    v = getvar(bd, var)
    return _interp1d_line(bd, v, point1, point2)
end

function _interp1d_line(bd::BatsrusIDLStructured{2, TV, TX, TW}, v, point1, point2) where {TV, TX, TW}
    xrange, yrange = get_range(bd)
    itp = linear_interp((xrange, yrange), parent(v))
    lx = point2[1] - point1[1]
    ly = point2[2] - point1[2]
    nx = lx ÷ step(xrange) |> Int
    ny = ly ÷ step(yrange) |> Int
    ns = floor(Int, √(nx^2 + ny^2))
    dx = lx / ns
    dy = ly / ns
    XQ = [point1[1] + i * dx for i in 0:ns]
    YQ = [point1[2] + i * dy for i in 0:ns]
    out = similar(XQ)
    itp(out, (XQ, YQ))

    return out
end

"""
    slice1d(bd, var, icut::Int=1, dir::Int=2)

Return view of variable `var` in `bd` along 1D slice. `icut` is the index along axis `dir`.
`dir == 1` means align with the 2nd (y) axis, `dir == 2` means align with the 1st (x) axis.
`var` can be an `AbstractString` or a `Symbol` for derived quantities.
"""
slice1d(bd, var::Union{AbstractString, Symbol}, icut::Int = 1, dir::Int = 2) =
    selectdim(getvar(bd, var), dir, icut)

"""
Return view of variable `var` in `bd`.
"""
@generated function getview(bd::BatsrusIDL{ndim}, var) where {ndim}
    colons = fill(:(:), ndim)
    return quote
        varIndex_ = findindex(bd, var)
        return @view bd.w[$(colons...), varIndex_]
    end
end

"""
Return value range of `var` in `bd`. `var` can be an `AbstractString` or a `Symbol`.
"""
get_var_range(bd::BatsrusIDL, var::Union{AbstractString, Symbol}) =
    getvar(bd, var) |> extrema

"""
Return mesh range of `bd`.
"""
@generated function get_range(bd::BatsrusIDLStructured{ndim}) where {ndim}
    return Expr(:tuple, [:(parent(val(dims(bd.x, $i)))) for i in 1:ndim]...)
end

"""
Squeeze singleton dimensions for an array `A`.
"""
function squeeze(A::AbstractArray)
    singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)

    return dropdims(A, dims = singleton_dims)
end

"""
    rotateTensorToVectorZ(tensor::AbstractMatrix, v::AbstractVector) -> SMatrix{3,3}

Rotate `tensor` with a rotation matrix that aligns the 3rd direction with `vector`, which is equivalent to change the basis from (i,j,k) to (i′,j′,k′) where k′ ∥ vector.
Reference: [Tensor rotation](https://math.stackexchange.com/questions/2303869/tensor-rotation)
"""
function rotateTensorToVectorZ(tensor::AbstractMatrix{T}, v) where {T}
    k = SVector{3, T}(0.0, 0.0, 1.0)
    axis = v × k |> normalize
    if axis[1] == axis[2] == 0
        return tensor
    else
        angle = acos(v ⋅ k / sqrt(v[1]^2 + v[2]^2 + v[3]^2))
        R = getRotationMatrix(axis, angle)
        return R * tensor * R'
    end
end

"""
    getRotationMatrix(axis::AbstractVector, angle::Real) --> SMatrix{3,3}

Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in radians.
Reference: [Rotation matrix from axis and angle](https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)

# Example

```julia
using LinearAlgebra
v = [-0.5, 1.0, 1.0]
v̂ = normalize(v)
θ = deg2rad(-74)
R = getRotationMatrix(v̂, θ)
```
"""
function getRotationMatrix(v::AbstractVector{<:AbstractFloat}, θ::Real)
    sinθ, cosθ = sincos(eltype(v)(θ))
    tmp = 1 - cosθ
    return m = @SMatrix [
        cosθ + v[1]^2 * tmp v[1] * v[2] * tmp - v[3] * sinθ v[1] * v[3] * tmp + v[2] * sinθ;
        v[1] * v[2] * tmp + v[3] * sinθ cosθ + v[2]^2 * tmp v[2] * v[3] * tmp - v[1] * sinθ;
        v[1] * v[3] * tmp - v[2] * sinθ v[3] * v[2] * tmp + v[1] * sinθ cosθ + v[3]^2 * tmp
    ]
end

"""
    generate_mock_amrex_data(output_dir::String; num_particles::Int=10, 
        real_component_names::Vector{String}=["u", "v"],
        particle_gen::Function)

Generate mock AMReX particle data for testing and benchmarking.

# Arguments

  - `output_dir::String`: Directory to save the mock data.
  - `num_particles::Int`: Number of particles to generate.
  - `real_component_names::Vector{String}`: Names of the extra real components (beyond x, y, z).
  - `particle_gen::Function`: A function `(i, n_reals) -> tuple` that takes an index `i` (1-based) and the total number of real components `n_reals`. It should return a tuple of `n_reals` Float64 values: `(x, y, z, comp1, comp2, ...)`.
"""
function generate_mock_amrex_data(
        output_dir::String;
        num_particles::Int = 10,
        real_component_names::Vector{String} = ["u", "v"],
        domain_min::Vector{Float64} = [0.0, 0.0, 0.0],
        domain_max::Vector{Float64} = [10.0, 10.0, 10.0],
        particle_gen::Function = (i, n_reals) -> (
            Float64(i), Float64(i), Float64(i), Float64(i * 10), Float64(i * 100),
        )
    )
    ptype = "particles"
    base_dir = joinpath(output_dir, ptype)
    mkpath(base_dir)

    n_extra = length(real_component_names)

    # Create Header
    header_path = joinpath(base_dir, "Header")
    open(header_path, "w") do f
        println(f, "Version_double")
        println(f, "3") # dim
        println(f, "$n_extra") # num_real_extra
        for name in real_component_names
            println(f, name)
        end
        println(f, "2") # num_int_extra (total int = 2 + 2 = 4)
        println(f, "id_1")
        println(f, "id_2")
        println(f, "0") # is_checkpoint (False)
        println(f, "$num_particles") # num_particles
        println(f, "$(num_particles + 1)") # max_next_id
        println(f, "0") # finest_level
        println(f, "1") # grids_per_level[0]
        println(f, "1 $num_particles 0") # grid info: which, count, where
    end

    # Create Level directory
    level_dir = joinpath(base_dir, "Level_0")
    mkpath(level_dir)

    # Create Particle_H
    particle_h_path = joinpath(level_dir, "Particle_H")
    open(particle_h_path, "w") do f
        # Use a single large box for simplicity
        # Convert physical range to index range (assuming dx=1 for simplicity if not specified)
        # Or just use arbitrary large indices since we only care about coverage
        # But wait, AMReXParticle uses domain_dimensions from Main Header to compute dx.
        # And Main Header uses domain indices.
        # Let's keep dx=1 effectively.
        dims = round.(Int, domain_max .- domain_min)
        # Ensure at least 1
        dims = max.(dims, 1)

        println(f, "(1 0") # num_boxes level
        # Box indices: lo hi (inclusive)
        # We just say the box covers the whole domain 0 to dims-1
        lo_str = "0,0,0"
        hi_str = join([d - 1 for d in dims], ",")
        println(f, "(($lo_str) ($hi_str) (0,0,0))")
    end

    # Create Data file
    data_fn = joinpath(level_dir, "DATA_00001")
    open(data_fn, "w") do f
        # Write particles
        # structure: num_particles * (3+n_extra) reals
        # x, y, z, ...
        n_reals = 3 + n_extra
        data = zeros(Float64, num_particles * n_reals)
        for i in 1:num_particles
            vals = particle_gen(i, n_reals)
            for j in 1:n_reals
                data[(i - 1) * n_reals + j] = vals[j]
            end
        end
        write(f, data)
    end

    # Create Main Header (for domain info)
    main_header_path = joinpath(output_dir, "Header")
    open(main_header_path, "w") do f
        println(f, "HyperCLaw-V1.1")
        println(f, "0") # num_fields
        println(f, "3") # dim
        dims = round.(Int, domain_max .- domain_min)
        dims = max.(dims, 1)
        println(f, "0.0") # time
        println(f, "0") # refine_ratio
        println(f, join(domain_min, " ")) # left_edge
        println(f, join(domain_max, " ")) # right_edge
        println(f, "0")

        lo_str = "0,0,0"
        hi_str = join([d - 1 for d in dims], ",")
        println(f, "(($lo_str) ($hi_str) (0,0,0))") # domain size
    end

    return
end

"""
    get_particle_field_aligned_transform(b_field, e_field=nothing)

Return a transformation function that converts particle data to a field-aligned coordinate system.

If only `b_field` is provided, the velocity components are decomposed into parallel and perpendicular to the magnetic field.
If `e_field` is also provided, an orthonormal basis (\$\\hat{\\mathbf{b}}\$, \$\\hat{\\mathbf{e}}\$, \$\\hat{\\mathbf{d}}\$) is created,
where \$\\hat{\\mathbf{d}} \\propto \\mathbf{B} \\times \\mathbf{E}\$ and \$\\hat{\\mathbf{e}} = \\hat{\\mathbf{d}} \\times \\hat{\\mathbf{b}}\$.

The returned function takes `(data, names)` and returns `(new_data, new_names)`.
"""
function get_particle_field_aligned_transform(b_field, e_field = nothing)
    bhat = SVector{3}(normalize(b_field))

    if isnothing(e_field)
        return (data, names) -> begin
            local idx_vx::Union{Nothing, Int}
            local idx_vy::Union{Nothing, Int}
            local idx_vz::Union{Nothing, Int}
            idx_vx = findfirst(x -> x in ("vx", "u", "ux", "velocity_x"), names)
            idx_vy = findfirst(x -> x in ("vy", "v", "uy", "velocity_y"), names)
            idx_vz = findfirst(x -> x in ("vz", "w", "uz", "velocity_z"), names)

            # Check for weights
            idx_w = findfirst(x -> x in ("weight", "w"), names)

            if isnothing(idx_vx) || isnothing(idx_vy) || isnothing(idx_vz)
                error("Velocity components not found in data.")
            end

            n_dims = size(data)
            n_names = length(names)

            if n_dims[1] != n_names
                error("Data dimensions $(n_dims) do not match component names length $(n_names) (expected (n_comp, n_part)).")
            end

            n_particles = n_dims[2]

            # output: v_para, v_perp, weight
            new_data = zeros(eltype(data), 3, n_particles)

            bx, by, bz = bhat
            @inbounds @simd for i in 1:n_particles
                vx = data[idx_vx, i]
                vy = data[idx_vy, i]
                vz = data[idx_vz, i]
                w = isnothing(idx_w) ? one(eltype(data)) : data[idx_w, i]

                v_para = vx * bx + vy * by + vz * bz

                vpx = vx - v_para * bx
                vpy = vy - v_para * by
                vpz = vz - v_para * bz
                v_perp = sqrt(vpx^2 + vpy^2 + vpz^2)

                new_data[1, i] = v_para
                new_data[2, i] = v_perp
                new_data[3, i] = w
            end

            return new_data, ["v_parallel", "v_perp", "weight"]
        end
    else # E and B provided
        # We want an orthonormal basis:
        # 1. b_hat
        # 2. e_perp_hat (along E perpendicular to B)
        # 3. exb_hat (along B x E)
        ehat_raw = SVector{3}(normalize(e_field))

        # Direction of B x E
        bxe = bhat × ehat_raw
        if norm(bxe) < 1.0e-9
            error("B and E are parallel, cannot define perpendicular directions uniquely.")
        end
        exb_hat = normalize(bxe)

        # Direction of E_perp (which is (B x E) x B)
        eperp_hat = normalize(exb_hat × bhat)

        (data, names) -> begin
            local idx_vx::Union{Nothing, Int}
            local idx_vy::Union{Nothing, Int}
            local idx_vz::Union{Nothing, Int}
            idx_vx = findfirst(x -> x in ("vx", "u", "ux", "velocity_x"), names)
            idx_vy = findfirst(x -> x in ("vy", "v", "uy", "velocity_y"), names)
            idx_vz = findfirst(x -> x in ("vz", "w", "uz", "velocity_z"), names)

            # Check for weights
            idx_w = findfirst(x -> x in ("weight", "w"), names)

            if isnothing(idx_vx) || isnothing(idx_vy) || isnothing(idx_vz)
                error("Velocity components not found in data.")
            end

            n_dims = size(data)
            n_names = length(names)

            if n_dims[1] != n_names
                error("Data dimensions $(n_dims) do not match component names length $(n_names) (expected (n_comp, n_part)).")
            end

            n_particles = n_dims[2]

            # output: v_B, v_E, v_BxE, weight
            new_data = zeros(eltype(data), 4, n_particles)

            bx, by, bz = bhat
            ex, ey, ez = eperp_hat
            dx, dy, dz = exb_hat

            @inbounds @simd for i in 1:n_particles
                vx = data[idx_vx, i]
                vy = data[idx_vy, i]
                vz = data[idx_vz, i]
                w = isnothing(idx_w) ? one(eltype(data)) : data[idx_w, i]

                # Project onto basis vectors
                v_b = vx * bx + vy * by + vz * bz
                v_e = vx * ex + vy * ey + vz * ez  # User defined v_E as along E_perp
                v_exb = vx * dx + vy * dy + vz * dz  # User defined v_BxE

                new_data[1, i] = v_b
                new_data[2, i] = v_e
                new_data[3, i] = v_exb
                new_data[4, i] = w
            end

            return new_data, ["v_B", "v_E", "v_BxE", "weight"]
        end
    end
end

function _triangulate_matplotlib end
