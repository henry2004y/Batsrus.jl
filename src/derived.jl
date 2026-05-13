# Derived quantities from raw output variables.
 
const FAC_J_PLANETARY = 10.0 / (4.0 * π * 6378.0)

"""
    get_magnitude2(bd::BatsrusIDL, var)

Calculate the magnitude square of vector `var`. See [`get_vectors`](@ref) for the options.
"""
function get_magnitude2(bd::BatsrusIDL{ndim}, var = :B) where {ndim}
    ivx, ivy, ivz = get_vectors_indices(bd, var)
    v = similar(bd.w, dims(bd.w)[1:ndim])
    w, v_raw = parent(bd.w), parent(v)

    @inbounds @simd for i in CartesianIndices(v_raw)
        v_raw[i] = w[i, ivx]^2 + w[i, ivy]^2 + w[i, ivz]^2
    end

    return v
end

"""
    get_magnitude(bd::BatsrusIDL, var)

Calculate the magnitude of vector `var`. See [`get_vectors`](@ref) for the options.
"""
function get_magnitude(bd::BatsrusIDL{ndim}, var = :B) where {ndim}
    ivx, ivy, ivz = get_vectors_indices(bd, var)
    v = similar(bd.w, dims(bd.w)[1:ndim])
    w, v_raw = parent(bd.w), parent(v)

    @inbounds @simd for i in CartesianIndices(v_raw)
        v_raw[i] = √(w[i, ivx]^2 + w[i, ivy]^2 + w[i, ivz]^2)
    end

    return v
end

"""
Return a tuple of indices for the vector of `var`. `var` can be `:B`, `:E`, `:U`, or any `:U`
followed by an index (e.g. `:U0` for species 0, `:U1` for species 1, etc.).
"""
function get_vectors_indices(bd::BatsrusIDL, var)
    if var === :B
        idx = findindex(bd, "bx")
    elseif var === :E
        idx = findindex(bd, "ex")
    elseif var === :U
        idx = findindex(bd, "ux")
    elseif var === :J
        idx = findindex(bd, "jx")
    else
        str = string(var)
        m = match(r"^U(\d+)$", str)
        if !isnothing(m)
            idx = findindex(bd, "uxs" * m[1])
        else
            throw(ArgumentError("Vector variable $var not supported"))
        end
    end
    return idx, idx + 1, idx + 2
end

"""
    get_vectors(bd::BatsrusIDL, var)

Return a tuple of vectors of `var`. `var` can be `:B`, `:E`, `:U`, or any `:U` followed by an
index (e.g. `:U0` for species 0, `:U1` for species 1, etc.).
"""
function get_vectors(bd::BatsrusIDL, var)
    ivx, ivy, ivz = get_vectors_indices(bd, var)
    return selectdim(bd.w, ndims(bd.w), ivx), selectdim(bd.w, ndims(bd.w), ivy),
        selectdim(bd.w, ndims(bd.w), ivz)
end

"""
    get_anisotropy(bd::BatsrusIDL, species=0)

Calculate the pressure anisotropy for `species`, indexing from 0. The default `method` is
based on the fact that the trace of the pressure tensor is a constant. The `rotation`
method is based on rotating the tensor.
"""
function get_anisotropy(bd::BatsrusIDL{2, TV}, species = 0; method = :simple) where {TV}
    ibx, iby, ibz = findindex(bd, "bx"), findindex(bd, "by"), findindex(bd, "bz")
    if species == 0
        ipxx = findindex(bd, "pxxs0")
    elseif species == 1
        ipxx = findindex(bd, "pxxs1")
    else
        ipxx = findindex(bd, "pxxs" * string(species))
    end
    ipyy, ipzz = ipxx + 1, ipxx + 2
    ipxy, ipxz, ipyz = ipxx + 3, ipxx + 4, ipxx + 5

    Paniso = similar(bd.w, dims(bd.w)[1:2])
    w, Paniso_raw = parent(bd.w), parent(Paniso)

    @inbounds for j in axes(w, 2), i in axes(w, 1)
        b̂ = normalize(SA[w[i, j, ibx], w[i, j, iby], w[i, j, ibz]])
        P = @SMatrix [
            w[i, j, ipxx] w[i, j, ipxy] w[i, j, ipxz];
            w[i, j, ipxy] w[i, j, ipyy] w[i, j, ipyz];
            w[i, j, ipxz] w[i, j, ipyz] w[i, j, ipzz]
        ]

        if method == :simple
            p_parallel = b̂' * P * b̂
            p_perp = (tr(P) - p_parallel) / 2
            Paniso_raw[i, j] = p_perp / p_parallel
        elseif method == :rotation
            Prot = rotateTensorToVectorZ(P, b̂)
            Paniso_raw[i, j] = (Prot[1, 1] + Prot[2, 2]) / (2 * Prot[3, 3])
        else
            error("Unknown method for get_anisotropy: $method. Use :simple or :rotation.")
        end
    end

    return Paniso
end

"""
    get_convection_E(bd::BatsrusIDL)

Return the convection electric field ``\\mathbf{E} = -\\mathbf{u}_i \\times \\mathbf{B}``.
"""
function get_convection_E(bd::BatsrusIDL)
    Bx, By, Bz = get_vectors(bd, :B)
    # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
    uix, uiy, uiz = get_vectors(bd, :U1)

    Econvx = similar(Bx)
    Econvy = similar(By)
    Econvz = similar(Bz)
    # -Ui × B
    @simd for i in eachindex(Econvx)
        Econvx[i] = -uiy[i] * Bz[i] + uiz[i] * By[i]
        Econvy[i] = -uiz[i] * Bx[i] + uix[i] * Bz[i]
        Econvz[i] = -uix[i] * By[i] + uiy[i] * Bx[i]
    end

    return Econvx, Econvy, Econvz
end

"""
    get_hall_E(bd::BatsrusIDL)

Return the Hall electric field ``\\mathbf{E}_H = (\\mathbf{u}_i - \\mathbf{u}_e) \\times \\mathbf{B}``.
"""
function get_hall_E(bd::BatsrusIDL)
    Bx, By, Bz = get_vectors(bd, :B)
    uex, uey, uez = get_vectors(bd, :U0)
    # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
    uix, uiy, uiz = get_vectors(bd, :U1)

    Ehallx = similar(Bx)
    Ehally = similar(By)
    Ehallz = similar(Bz)
    # (Ui - Ue) × B
    for i in eachindex(Ehallx)
        Ehallx[i] = (uiy[i] - uey[i]) * Bz[i] - (uiz[i] - uez[i]) * By[i]
        Ehally[i] = (uiz[i] - uez[i]) * Bx[i] - (uix[i] - uex[i]) * Bz[i]
        Ehallz[i] = (uix[i] - uex[i]) * By[i] - (uiy[i] - uey[i]) * Bx[i]
    end

    return Ehallx, Ehally, Ehallz
end

"""
    fill_vector_from_scalars(bd::BatsrusIDL, var)

Construct vector of `var` from its scalar components. Alternatively, check
[`get_vectors`](@ref) for returning vector components as separate arrays.
"""
function fill_vector_from_scalars(bd::BatsrusIDL, var)
    vt = get_vectors(bd, var)
    Rpost = CartesianIndices(size(bd.x)[1:(end - 1)])
    return v = @inbounds [vt[iv][i] for iv in 1:3, i in Rpost]
end

"""
    get_current_density(bd::BatsrusIDL)

Calculate the current density ``\\mathbf{J} = \\nabla \\times \\mathbf{B} / \\mu_0`` from the curl of
the magnetic field. For PLANETARY units, the result is in ``\\mu A/m^2``. Currently only
supports structured grids.
"""
function get_current_density(bd::BatsrusIDLStructured{1, TV}) where {TV}
    ivx, ivy, ivz = get_vectors_indices(bd, :B)
    w = parent(bd.w)
    xrange = get_range(bd)[1]
    dx = TV(step(xrange))
    nx = size(w, 1)

    jx = zeros(TV, nx)
    jy = zeros(TV, nx)
    jz = zeros(TV, nx)

    # J = ∇ × B. In 1D, ∂/∂y = ∂/∂z = 0.
    # Jx = ∂Bz/∂y - ∂By/∂z = 0
    # Jy = ∂Bx/∂z - ∂Bz/∂x = -∂Bz/∂x
    # Jz = ∂By/∂x - ∂Bx/∂y = ∂By/∂x

    @inbounds for ix in 1:nx
        if ix == 1
            dbz_dx = (w[2, ivz] - w[1, ivz]) / dx
            dby_dx = (w[2, ivy] - w[1, ivy]) / dx
        elseif ix == nx
            dbz_dx = (w[nx, ivz] - w[nx - 1, ivz]) / dx
            dby_dx = (w[nx, ivy] - w[nx - 1, ivy]) / dx
        else
            dbz_dx = (w[ix + 1, ivz] - w[ix - 1, ivz]) / (2dx)
            dby_dx = (w[ix + 1, ivy] - w[ix - 1, ivy]) / (2dx)
        end
        jy[ix] = -dbz_dx
        jz[ix] = dby_dx
    end

    if bd.head.headline == "PLANETARY"
        fac = TV(FAC_J_PLANETARY)
        jy .*= fac
        jz .*= fac
    end

    # Wrap in DimArrays
    dims_ = dims(bd.w)[1:1]
    Jx = DimArray(jx, dims_)
    Jy = DimArray(jy, dims_)
    Jz = DimArray(jz, dims_)

    return Jx, Jy, Jz
end

function get_current_density(bd::BatsrusIDLStructured{2, TV}) where {TV}
    ivx, ivy, ivz = get_vectors_indices(bd, :B)
    w = parent(bd.w)
    xrange, yrange = get_range(bd)
    dx = TV(step(xrange))
    dy = TV(step(yrange))
    nx, ny = size(w, 1), size(w, 2)

    jx = zeros(TV, nx, ny)
    jy = zeros(TV, nx, ny)
    jz = zeros(TV, nx, ny)

    # Jx = ∂Bz/∂y - ∂By/∂z = ∂Bz/∂y
    # Jy = ∂Bx/∂z - ∂Bz/∂x = -∂Bz/∂x
    # Jz = ∂By/∂x - ∂Bx/∂y

    @inbounds for iy in 1:ny, ix in 1:nx
        # ∂Bz/∂y
        if iy == 1
            dbz_dy = (w[ix, 2, ivz] - w[ix, 1, ivz]) / dy
        elseif iy == ny
            dbz_dy = (w[ix, ny, ivz] - w[ix, ny - 1, ivz]) / dy
        else
            dbz_dy = (w[ix, iy + 1, ivz] - w[ix, iy - 1, ivz]) / (2dy)
        end

        # ∂Bz/∂x
        if ix == 1
            dbz_dx = (w[2, iy, ivz] - w[1, iy, ivz]) / dx
        elseif ix == nx
            dbz_dx = (w[nx, iy, ivz] - w[nx - 1, iy, ivz]) / dx
        else
            dbz_dx = (w[ix + 1, iy, ivz] - w[ix - 1, iy, ivz]) / (2dx)
        end

        # ∂By/∂x
        if ix == 1
            dby_dx = (w[2, iy, ivy] - w[1, iy, ivy]) / dx
        elseif ix == nx
            dby_dx = (w[nx, iy, ivy] - w[nx - 1, iy, ivy]) / dx
        else
            dby_dx = (w[ix + 1, iy, ivy] - w[ix - 1, iy, ivy]) / (2dx)
        end

        # ∂Bx/∂y
        if iy == 1
            dbx_dy = (w[ix, 2, ivx] - w[ix, 1, ivx]) / dy
        elseif iy == ny
            dbx_dy = (w[ix, ny, ivx] - w[ix, ny - 1, ivx]) / dy
        else
            dbx_dy = (w[ix, iy + 1, ivx] - w[ix, iy - 1, ivx]) / (2dy)
        end

        jx[ix, iy] = dbz_dy
        jy[ix, iy] = -dbz_dx
        jz[ix, iy] = dby_dx - dbx_dy
    end

    if bd.head.headline == "PLANETARY"
        fac = TV(FAC_J_PLANETARY)
        jx .*= fac
        jy .*= fac
        jz .*= fac
    end

    dims_ = dims(bd.w)[1:2]
    Jx = DimArray(jx, dims_)
    Jy = DimArray(jy, dims_)
    Jz = DimArray(jz, dims_)

    return Jx, Jy, Jz
end

function get_current_density(bd::BatsrusIDLStructured{3, TV}) where {TV}
    ivx, ivy, ivz = get_vectors_indices(bd, :B)
    w = parent(bd.w)
    xrange, yrange, zrange = get_range(bd)
    dx = TV(step(xrange))
    dy = TV(step(yrange))
    dz = TV(step(zrange))
    nx, ny, nz = size(w, 1), size(w, 2), size(w, 3)

    jx = zeros(TV, nx, ny, nz)
    jy = zeros(TV, nx, ny, nz)
    jz = zeros(TV, nx, ny, nz)

    @inbounds for iz in 1:nz, iy in 1:ny, ix in 1:nx
        # ∂Bz/∂y
        if iy == 1
            dbz_dy = (w[ix, 2, iz, ivz] - w[ix, 1, iz, ivz]) / dy
        elseif iy == ny
            dbz_dy = (w[ix, ny, iz, ivz] - w[ix, ny - 1, iz, ivz]) / dy
        else
            dbz_dy = (w[ix, iy + 1, iz, ivz] - w[ix, iy - 1, iz, ivz]) / (2dy)
        end

        # ∂By/∂z
        if iz == 1
            dby_dz = (w[ix, iy, 2, ivy] - w[ix, iy, 1, ivy]) / dz
        elseif iz == nz
            dby_dz = (w[ix, iy, nz, ivy] - w[ix, iy, nz - 1, ivy]) / dz
        else
            dby_dz = (w[ix, iy, iz + 1, ivy] - w[ix, iy, iz - 1, ivy]) / (2dz)
        end

        # ∂Bx/∂z
        if iz == 1
            dbx_dz = (w[ix, iy, 2, ivx] - w[ix, iy, 1, ivx]) / dz
        elseif iz == nz
            dbx_dz = (w[ix, iy, nz, ivx] - w[ix, iy, nz - 1, ivx]) / dz
        else
            dbx_dz = (w[ix, iy, iz + 1, ivx] - w[ix, iy, iz - 1, ivx]) / (2dz)
        end

        # ∂Bz/∂x
        if ix == 1
            dbz_dx = (w[2, iy, iz, ivz] - w[1, iy, iz, ivz]) / dx
        elseif ix == nx
            dbz_dx = (w[nx, iy, iz, ivz] - w[nx - 1, iy, iz, ivz]) / dx
        else
            dbz_dx = (w[ix + 1, iy, iz, ivz] - w[ix - 1, iy, iz, ivz]) / (2dx)
        end

        # ∂By/∂x
        if ix == 1
            dby_dx = (w[2, iy, iz, ivy] - w[1, iy, iz, ivy]) / dx
        elseif ix == nx
            dby_dx = (w[nx, iy, iz, ivy] - w[nx - 1, iy, iz, ivy]) / dx
        else
            dby_dx = (w[ix + 1, iy, iz, ivy] - w[ix - 1, iy, iz, ivy]) / (2dx)
        end

        # ∂Bx/∂y
        if iy == 1
            dbx_dy = (w[ix, 2, iz, ivx] - w[ix, 1, iz, ivx]) / dy
        elseif iy == ny
            dbx_dy = (w[ix, ny, iz, ivx] - w[ix, ny - 1, iz, ivx]) / dy
        else
            dbx_dy = (w[ix, iy + 1, iz, ivx] - w[ix, iy - 1, iz, ivx]) / (2dy)
        end

        jx[ix, iy, iz] = dbz_dy - dby_dz
        jy[ix, iy, iz] = dbx_dz - dbz_dx
        jz[ix, iy, iz] = dby_dx - dbx_dy
    end

    if bd.head.headline == "PLANETARY"
        fac = TV(FAC_J_PLANETARY)
        jx .*= fac
        jy .*= fac
        jz .*= fac
    end

    dims_ = dims(bd.w)[1:3]
    Jx = DimArray(jx, dims_)
    Jy = DimArray(jy, dims_)
    Jz = DimArray(jz, dims_)

    return Jx, Jy, Jz
end

function get_current_density(bd::BatsrusIDLUnstructured)
    error("get_current_density is not yet supported for unstructured grids.")
end

# --- Derived scalar quantities ---
@inline _getvar(bd::BatsrusIDL, ::Val{:b}) = get_magnitude(bd, :B)
@inline _getvar(bd::BatsrusIDL, ::Val{:b2}) = get_magnitude2(bd, :B)
@inline _getvar(bd::BatsrusIDL, ::Val{:e}) = get_magnitude(bd, :E)
@inline _getvar(bd::BatsrusIDL, ::Val{:u}) = get_magnitude(bd, :U)
@inline _getvar(bd::BatsrusIDL{2}, ::Val{:anisotropy0}) = get_anisotropy(bd, 0)
@inline _getvar(bd::BatsrusIDL{2}, ::Val{:anisotropy1}) = get_anisotropy(bd, 1)

@inline function _getvar(bd::BatsrusIDL, ::Val{:jx})
    Jx, _, _ = get_current_density(bd)
    return Jx
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:jy})
    _, Jy, _ = get_current_density(bd)
    return Jy
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:jz})
    _, _, Jz = get_current_density(bd)
    return Jz
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:j})
    Jx, Jy, Jz = get_current_density(bd)
    return sqrt.(Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2)
end
