# Derived quantities from raw output variables.

# Earth radius in km used for current density scaling in PLANETARY units.
const EARTH_RADIUS_KM = 6378.0
# Factor to convert curl(B) to current density in μA/m²: 10 / (4π * Re)
const FAC_J_PLANETARY = 10.0 / (4.0 * π * EARTH_RADIUS_KM)


"""
    get_magnitude2(bd::BatsrusIDL, var)

Calculate the magnitude square of vector `var`. See [`get_vectors`](@ref) for the options.
"""
function get_magnitude2(bd::BatsrusIDL, var = :B)
    vx, vy, vz = get_vectors(bd, var)
    return vx .^ 2 .+ vy .^ 2 .+ vz .^ 2
end

"""
    get_magnitude(bd::BatsrusIDL, var = :B)

Calculate the magnitude of vector `var`. See [`get_vectors`](@ref) for the options.
"""
@inline function get_magnitude(bd::BatsrusIDL{ndim, TV}, var = :B) where {ndim, TV}
    _get_magnitude(bd, Val(var))
end

@inline function _get_magnitude(bd::BatsrusIDLStructured{ndim, TV, TX, TW}, ::Val{var}) where {ndim, TV, TX, TW, var}
    indices = get_vectors_indices(bd, Val(var))
    w = parent(bd.w)
    
    ivx, ivy, ivz = indices
    d = dims(bd.w)
    
    if ndim == 2
        nx, ny = size(bd.w, 1), size(bd.w, 2)
        n_space = nx * ny
        res = similar(w, nx, ny)
        @inbounds for i in 1:n_space
            res[i] = sqrt(w[i + (ivx-1)*n_space]^2 + 
                          w[i + (ivy-1)*n_space]^2 + 
                          w[i + (ivz-1)*n_space]^2)
        end
        return rebuild(bd.w, res, (d[1], d[2]))
    elseif ndim == 3
        nx, ny, nz = size(bd.w, 1), size(bd.w, 2), size(bd.w, 3)
        n_space = nx * ny * nz
        res = similar(w, nx, ny, nz)
        @inbounds for i in 1:n_space
            res[i] = sqrt(w[i + (ivx-1)*n_space]^2 + 
                          w[i + (ivy-1)*n_space]^2 + 
                          w[i + (ivz-1)*n_space]^2)
        end
        return rebuild(bd.w, res, (d[1], d[2], d[3]))
    else
        sz = ntuple(i -> size(bd.w, i), Val(ndim))
        n_space = 1
        for i in 1:ndim; n_space *= sz[i]; end
        res = similar(w, sz)
        @inbounds for i in 1:n_space
            res[i] = sqrt(w[i + (ivx-1)*n_space]^2 + 
                          w[i + (ivy-1)*n_space]^2 + 
                          w[i + (ivz-1)*n_space]^2)
        end
        return rebuild(bd.w, res, ntuple(i -> d[i], Val(ndim)))
    end
end

@inline function _get_magnitude(bd::BatsrusIDLUnstructured{ndim, TV, TX, TW}, ::Val{var}) where {ndim, TV, TX, TW, var}
    indices = get_vectors_indices(bd, Val(var))
    w = parent(bd.w)
    n_cells = size(w, 1)
    
    res = similar(w, n_cells)
    ivx, ivy, ivz = indices
    
    @inbounds for i in 1:n_cells
        res[i] = sqrt(w[i + (ivx-1)*n_cells]^2 + 
                      w[i + (ivy-1)*n_cells]^2 + 
                      w[i + (ivz-1)*n_cells]^2)
    end
    
    d = dims(bd.w)
    return DimArray(res, (d[1],))
end

@inline function _get_pressure_tensor_indices(bd::BatsrusIDL, species)
    if species == 0
        idx = findindex(bd, "pxxs0")
        return (idx, idx + 1, idx + 2, idx + 3, idx + 4, idx + 5)
    elseif species == 1
        idx = findindex(bd, "pxxs1")
        return (idx, idx + 1, idx + 2, idx + 3, idx + 4, idx + 5)
    elseif species == 2
        idx = findindex(bd, "pxxs2")
        return (idx, idx + 1, idx + 2, idx + 3, idx + 4, idx + 5)
    else
        idx = findindex(bd, "pxxs" * string(species))
        return (idx, idx + 1, idx + 2, idx + 3, idx + 4, idx + 5)
    end
end

"""
    get_anisotropy(bd::BatsrusIDL, species=0; method=:projection)

Calculate the pressure anisotropy for `species`. Two methods are supported:
- `:projection`: direct projection of the pressure tensor onto the magnetic field (default).
- `:rotation`: rotating the pressure tensor to a field-aligned coordinate system.
"""
@inline function get_anisotropy(bd::BatsrusIDL{ndim, TV}, species = 0; method = :projection) where {ndim, TV}
    _get_anisotropy(bd, species, method)
end

@inline function _get_anisotropy(bd::BatsrusIDLStructured{ndim, TV, TX, TW}, species, method) where {ndim, TV, TX, TW}
    iv = get_vectors_indices(bd, :B)
    ip = _get_pressure_tensor_indices(bd, species)
    
    w = parent(bd.w)
    d = dims(bd.w)

    if ndim == 2
        nx, ny = size(bd.w, 1), size(bd.w, 2)
        n_space = nx * ny
        res = similar(w, nx, ny)
        
        if method === :projection
            ivx, ivy, ivz = iv
            ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip
            @inbounds for i in 1:n_space
                Bx, By, Bz = w[i + (ivx-1)*n_space], w[i + (ivy-1)*n_space], w[i + (ivz-1)*n_space]
                B2 = Bx^2 + By^2 + Bz^2
                pxx, pyy, pzz = w[i + (ipxx-1)*n_space], w[i + (ipyy-1)*n_space], w[i + (ipzz-1)*n_space]
                pxy, pxz, pyz = w[i + (ipxy-1)*n_space], w[i + (ipxz-1)*n_space], w[i + (ipyz-1)*n_space]
                p_parallel = (pxx*Bx^2 + pyy*By^2 + pzz*Bz^2 + 2*pxy*Bx*By + 2*pxz*Bx*Bz + 2*pyz*By*Bz) / B2
                res[i] = (pxx + pyy + pzz - p_parallel) / (2 * p_parallel)
            end
        elseif method === :rotation
            ivx, ivy, ivz = iv
            ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip
            @inbounds for i in 1:n_space
                Bx, By, Bz = w[i + (ivx-1)*n_space], w[i + (ivy-1)*n_space], w[i + (ivz-1)*n_space]
                pxx, pyy, pzz = w[i + (ipxx-1)*n_space], w[i + (ipyy-1)*n_space], w[i + (ipzz-1)*n_space]
                pxy, pxz, pyz = w[i + (ipxy-1)*n_space], w[i + (ipxz-1)*n_space], w[i + (ipyz-1)*n_space]
                P = @SMatrix [pxx pxy pxz; pxy pyy pyz; pxz pyz pzz]
                v = @SVector [Bx, By, Bz]
                Pr = rotateTensorToVectorZ(P, v)
                res[i] = (Pr[1, 1] + Pr[2, 2]) / (2 * Pr[3, 3])
            end
        else
            error("Unknown method $method")
        end
        return rebuild(bd.w, res, (d[1], d[2]))
    else
        # Fallback for other dimensions
        sz = ntuple(i -> size(bd.w, i), Val(ndim))
        n_space = 1
        for i in 1:ndim; n_space *= sz[i]; end
        res = similar(w, sz)
        # ...
        return rebuild(bd.w, res, ntuple(i -> d[i], Val(ndim)))
    end
end

@inline function _get_anisotropy(bd::BatsrusIDLUnstructured{ndim, TV, TX, TW}, species, method) where {ndim, TV, TX, TW}
    iv = get_vectors_indices(bd, :B)
    ip = _get_pressure_tensor_indices(bd, species)
    
    w = parent(bd.w)
    n_cells = size(w, 1)
    
    res = similar(w, n_cells)
    
    if method === :projection
        ivx, ivy, ivz = iv
        ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip
        @inbounds for i in 1:n_cells
            Bx = w[i + (ivx-1)*n_cells]
            By = w[i + (ivy-1)*n_cells]
            Bz = w[i + (ivz-1)*n_cells]
            B2 = Bx^2 + By^2 + Bz^2
            
            pxx = w[i + (ipxx-1)*n_cells]
            pyy = w[i + (ipyy-1)*n_cells]
            pzz = w[i + (ipzz-1)*n_cells]
            pxy = w[i + (ipxy-1)*n_cells]
            pxz = w[i + (ipxz-1)*n_cells]
            pyz = w[i + (ipyz-1)*n_cells]
            
            p_parallel = (pxx*Bx^2 + pyy*By^2 + pzz*Bz^2 + 
                          2*pxy*Bx*By + 2*pxz*Bx*Bz + 2*pyz*By*Bz) / B2
            p_perp = (pxx + pyy + pzz - p_parallel) / 2
            res[i] = p_perp / p_parallel
        end
    elseif method === :rotation
        ivx, ivy, ivz = iv
        ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip
        @inbounds for i in 1:n_cells
            Bx = w[i + (ivx-1)*n_cells]
            By = w[i + (ivy-1)*n_cells]
            Bz = w[i + (ivz-1)*n_cells]
            
            pxx = w[i + (ipxx-1)*n_cells]
            pyy = w[i + (ipyy-1)*n_cells]
            pzz = w[i + (ipzz-1)*n_cells]
            pxy = w[i + (ipxy-1)*n_cells]
            pxz = w[i + (ipxz-1)*n_cells]
            pyz = w[i + (ipyz-1)*n_cells]

            P = @SMatrix [pxx pxy pxz; pxy pyy pyz; pxz pyz pzz]
            v = @SVector [Bx, By, Bz]
            Pr = rotateTensorToVectorZ(P, v)
            p_parallel = Pr[3, 3]
            p_perp = (Pr[1, 1] + Pr[2, 2]) / 2
            res[i] = p_perp / p_parallel
        end
    else
        error("Unknown method $method")
    end

    d = dims(bd.w)
    return DimArray(res, (d[1],))
end

"""
    get_convection_E(bd::BatsrusIDL)

Return the convection electric field ``\\mathbf{E} = -\\mathbf{u}_i \\times \\mathbf{B}``.
"""
function get_convection_E(bd::BatsrusIDL)
    Bx, By, Bz = get_vectors(bd, :B)
    # Fallback to species 1 if 'ux' is not found
    Utype = _has_var(bd, "ux") ? :U : :U1
    uix, uiy, uiz = get_vectors(bd, Utype)

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
    return @inbounds [vt[iv][i] for iv in 1:3, i in Rpost]
end

"""
    get_current_density(bd::BatsrusIDL)

Calculate the current density ``\\mathbf{J} = \\nabla \\times \\mathbf{B} / \\mu_0`` from the curl of
the magnetic field. For PLANETARY units, the result is in ``\\mu A/m^2``. Currently only
supports structured grids.
"""
function get_current_density(bd::BatsrusIDLStructured{1, TV}) where {TV}
    if _has_var(bd, "jx") && _has_var(bd, "jy") && _has_var(bd, "jz")
        return getvar(bd, "jx"), getvar(bd, "jy"), getvar(bd, "jz")
    end
    ivz = findindex(bd, "bz")
    ivy = findindex(bd, "by")
    w = parent(bd.w)
    xrange = get_range(bd)[1]
    dx = TV(step(xrange))
    nx = size(w, 1)

    jx = zeros(TV, nx)
    jy = zeros(TV, nx)
    jz = zeros(TV, nx)

    @inbounds for ix in 1:nx
        jy[ix] = -_diff1(w, ix, nx, dx, ivz)
        jz[ix] = _diff1(w, ix, nx, dx, ivy)
    end

    _apply_j_scaling!(jx, bd)
    _apply_j_scaling!(jy, bd)
    _apply_j_scaling!(jz, bd)

    dims_ = dims(bd.w)[1:1]
    return DimArray(jx, dims_), DimArray(jy, dims_), DimArray(jz, dims_)
end

function get_current_density(bd::BatsrusIDLStructured{2, TV}) where {TV}
    if _has_var(bd, "jx") && _has_var(bd, "jy") && _has_var(bd, "jz")
        return getvar(bd, "jx"), getvar(bd, "jy"), getvar(bd, "jz")
    end
    ivx, ivy, ivz = get_vectors_indices(bd, :B)
    w = parent(bd.w)
    xrange, yrange = get_range(bd)
    dx, dy = TV(step(xrange)), TV(step(yrange))
    nx, ny = size(w, 1), size(w, 2)

    jx = zeros(TV, nx, ny)
    jy = zeros(TV, nx, ny)
    jz = zeros(TV, nx, ny)

    @inbounds for iy in 1:ny, ix in 1:nx
        dbz_dy = _diff2_y(w, ix, iy, ny, dy, ivz)
        dbz_dx = _diff2_x(w, ix, iy, nx, dx, ivz)
        dby_dx = _diff2_x(w, ix, iy, nx, dx, ivy)
        dbx_dy = _diff2_y(w, ix, iy, ny, dy, ivx)

        jx[ix, iy] = dbz_dy
        jy[ix, iy] = -dbz_dx
        jz[ix, iy] = dby_dx - dbx_dy
    end

    _apply_j_scaling!(jx, bd)
    _apply_j_scaling!(jy, bd)
    _apply_j_scaling!(jz, bd)

    dims_ = dims(bd.w)[1:2]
    return DimArray(jx, dims_), DimArray(jy, dims_), DimArray(jz, dims_)
end

function get_current_density(bd::BatsrusIDLStructured{3, TV}) where {TV}
    if _has_var(bd, "jx") && _has_var(bd, "jy") && _has_var(bd, "jz")
        return getvar(bd, "jx"), getvar(bd, "jy"), getvar(bd, "jz")
    end
    ivx, ivy, ivz = get_vectors_indices(bd, :B)
    w = parent(bd.w)
    xrange, yrange, zrange = get_range(bd)
    dx, dy, dz = TV(step(xrange)), TV(step(yrange)), TV(step(zrange))
    nx, ny, nz = size(w, 1), size(w, 2), size(w, 3)

    jx = zeros(TV, nx, ny, nz)
    jy = zeros(TV, nx, ny, nz)
    jz = zeros(TV, nx, ny, nz)

    @inbounds for iz in 1:nz, iy in 1:ny, ix in 1:nx
        dbz_dy = _diff3_y(w, ix, iy, iz, ny, dy, ivz)
        dby_dz = _diff3_z(w, ix, iy, iz, nz, dz, ivy)
        dbx_dz = _diff3_z(w, ix, iy, iz, nz, dz, ivx)
        dbz_dx = _diff3_x(w, ix, iy, iz, nx, dx, ivz)
        dby_dx = _diff3_x(w, ix, iy, iz, nx, dx, ivy)
        dbx_dy = _diff3_y(w, ix, iy, iz, ny, dy, ivx)

        jx[ix, iy, iz] = dbz_dy - dby_dz
        jy[ix, iy, iz] = dbx_dz - dbz_dx
        jz[ix, iy, iz] = dby_dx - dbx_dy
    end

    _apply_j_scaling!(jx, bd)
    _apply_j_scaling!(jy, bd)
    _apply_j_scaling!(jz, bd)

    dims_ = dims(bd.w)[1:3]
    return DimArray(jx, dims_), DimArray(jy, dims_), DimArray(jz, dims_)
end

function get_current_density(bd::BatsrusIDLUnstructured)
    error("get_current_density is not yet supported for unstructured grids.")
end

# --- Derived scalar quantities ---

@inline function _getvar(bd::BatsrusIDL, ::Val{:jx})
    if _has_var(bd, "jx")
        return selectdim(bd.w, ndims(bd.w), findindex(bd, "jx"))
    end
    return _compute_jx(bd)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:jy})
    if _has_var(bd, "jy")
        return selectdim(bd.w, ndims(bd.w), findindex(bd, "jy"))
    end
    return _compute_jy(bd)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:jz})
    if _has_var(bd, "jz")
        return selectdim(bd.w, ndims(bd.w), findindex(bd, "jz"))
    end
    return _compute_jz(bd)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:j})
    if _has_var(bd, "jx") && _has_var(bd, "jy") && _has_var(bd, "jz")
        Jx = selectdim(bd.w, ndims(bd.w), findindex(bd, "jx"))
        Jy = selectdim(bd.w, ndims(bd.w), findindex(bd, "jy"))
        Jz = selectdim(bd.w, ndims(bd.w), findindex(bd, "jz"))
        return sqrt.(Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2)
    end
    Jx, Jy, Jz = get_current_density(bd)
    return sqrt.(Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2)
end

@inline _getvar(bd::BatsrusIDL, ::Val{:b}) = get_magnitude(bd, :B)
@inline _getvar(bd::BatsrusIDL, ::Val{:b2}) = get_magnitude2(bd, :B)
@inline _getvar(bd::BatsrusIDL, ::Val{:e}) = get_magnitude(bd, :E)
@inline _getvar(bd::BatsrusIDL, ::Val{:u}) = get_magnitude(bd, :U)
@inline _getvar(bd::BatsrusIDL{2}, ::Val{:anisotropy0}) = get_anisotropy(bd, 0)
@inline _getvar(bd::BatsrusIDL{2}, ::Val{:anisotropy1}) = get_anisotropy(bd, 1)

# --- Internal computation helpers ---

@inline _compute_jx(bd::BatsrusIDL) = error("jx computation not supported for this grid type.")
@inline _compute_jy(bd::BatsrusIDL) = error("jy computation not supported for this grid type.")
@inline _compute_jz(bd::BatsrusIDL) = error("jz computation not supported for this grid type.")

@inline function _compute_jx(bd::BatsrusIDLStructured{1, TV}) where {TV}
    return DimArray(zeros(TV, size(bd.w, 1)), dims(bd.w)[1:1])
end

@inline function _compute_jy(bd::BatsrusIDLStructured{1, TV}) where {TV}
    ivz = findindex(bd, "bz")
    w = parent(bd.w)
    xrange = get_range(bd)[1]
    dx, nx = TV(step(xrange)), size(w, 1)
    jy = [_diff1(w, ix, nx, dx, ivz) for ix in 1:nx]
    jy .*= -1 # Jy = -∂Bz/∂x
    return DimArray(_apply_j_scaling!(jy, bd), dims(bd.w)[1:1])
end

@inline function _compute_jz(bd::BatsrusIDLStructured{1, TV}) where {TV}
    ivy = findindex(bd, "by")
    w = parent(bd.w)
    xrange = get_range(bd)[1]
    dx, nx = TV(step(xrange)), size(w, 1)
    jz = [_diff1(w, ix, nx, dx, ivy) for ix in 1:nx]
    return DimArray(_apply_j_scaling!(jz, bd), dims(bd.w)[1:1])
end

@inline function _compute_jx(bd::BatsrusIDLStructured{2, TV}) where {TV}
    ivz = findindex(bd, "bz")
    w = parent(bd.w)
    yrange = get_range(bd)[2]
    dy, nx, ny = TV(step(yrange)), size(w, 1), size(w, 2)
    jx = [ _diff2_y(w, ix, iy, ny, dy, ivz) for ix in 1:nx, iy in 1:ny ]
    return DimArray(_apply_j_scaling!(jx, bd), dims(bd.w)[1:2])
end

@inline function _compute_jy(bd::BatsrusIDLStructured{2, TV}) where {TV}
    ivz = findindex(bd, "bz")
    w = parent(bd.w)
    xrange = get_range(bd)[1]
    dx, nx, ny = TV(step(xrange)), size(w, 1), size(w, 2)
    jy = [ -_diff2_x(w, ix, iy, nx, dx, ivz) for ix in 1:nx, iy in 1:ny ]
    return DimArray(_apply_j_scaling!(jy, bd), dims(bd.w)[1:2])
end

@inline function _compute_jz(bd::BatsrusIDLStructured{2, TV}) where {TV}
    ivx, ivy = findindex(bd, "bx"), findindex(bd, "by")
    w = parent(bd.w)
    xrange, yrange = get_range(bd)
    dx, dy = TV(step(xrange)), TV(step(yrange))
    nx, ny = size(w, 1), size(w, 2)
    jz = [
        _diff2_x(w, ix, iy, nx, dx, ivy) - _diff2_y(w, ix, iy, ny, dy, ivx)
            for ix in 1:nx, iy in 1:ny
    ]
    return DimArray(_apply_j_scaling!(jz, bd), dims(bd.w)[1:2])
end

@inline function _compute_jx(bd::BatsrusIDLStructured{3, TV}) where {TV}
    ivy, ivz = findindex(bd, "by"), findindex(bd, "bz")
    w = parent(bd.w)
    _, yrange, zrange = get_range(bd)
    dy, dz = TV(step(yrange)), TV(step(zrange))
    nx, ny, nz = size(w, 1), size(w, 2), size(w, 3)
    jx = [
        _diff3_y(w, ix, iy, iz, ny, dy, ivz) - _diff3_z(w, ix, iy, iz, nz, dz, ivy)
            for ix in 1:nx, iy in 1:ny, iz in 1:nz
    ]
    return DimArray(_apply_j_scaling!(jx, bd), dims(bd.w)[1:3])
end

@inline function _compute_jy(bd::BatsrusIDLStructured{3, TV}) where {TV}
    ivx, ivz = findindex(bd, "bx"), findindex(bd, "bz")
    w = parent(bd.w)
    xrange, _, zrange = get_range(bd)
    dx, dz = TV(step(xrange)), TV(step(zrange))
    nx, ny, nz = size(w, 1), size(w, 2), size(w, 3)
    jy = [
        _diff3_z(w, ix, iy, iz, nz, dz, ivx) - _diff3_x(w, ix, iy, iz, nx, dx, ivz)
            for ix in 1:nx, iy in 1:ny, iz in 1:nz
    ]
    return DimArray(_apply_j_scaling!(jy, bd), dims(bd.w)[1:3])
end

@inline function _compute_jz(bd::BatsrusIDLStructured{3, TV}) where {TV}
    ivx, ivy = findindex(bd, "bx"), findindex(bd, "by")
    w = parent(bd.w)
    xrange, yrange = get_range(bd)
    dx, dy = TV(step(xrange)), TV(step(yrange))
    nx, ny, nz = size(w, 1), size(w, 2), size(w, 3)
    jz = [
        _diff3_x(w, ix, iy, iz, nx, dx, ivy) - _diff3_y(w, ix, iy, iz, ny, dy, ivx)
            for ix in 1:nx, iy in 1:ny, iz in 1:nz
    ]
    return DimArray(_apply_j_scaling!(jz, bd), dims(bd.w)[1:3])
end

# --- Internal finite difference helpers ---

@inline function _apply_j_scaling!(j, bd::BatsrusIDL)
    if bd.head.headline == "PLANETARY"
        j .*= eltype(j)(FAC_J_PLANETARY)
    end
    return j
end

@inline function _diff1(w, i, n, h, iv)
    if i == 1
        return (w[2, iv] - w[1, iv]) / h
    elseif i == n
        return (w[n, iv] - w[n - 1, iv]) / h
    else
        return (w[i + 1, iv] - w[i - 1, iv]) / (2h)
    end
end

@inline function _diff2_x(w, ix, iy, nx, dx, iv)
    if ix == 1
        return (w[2, iy, iv] - w[1, iy, iv]) / dx
    elseif ix == nx
        return (w[nx, iy, iv] - w[nx - 1, iy, iv]) / dx
    else
        return (w[ix + 1, iy, iv] - w[ix - 1, iy, iv]) / (2dx)
    end
end

@inline function _diff2_y(w, ix, iy, ny, dy, iv)
    if iy == 1
        return (w[ix, 2, iv] - w[ix, 1, iv]) / dy
    elseif iy == ny
        return (w[ix, ny, iv] - w[ix, ny - 1, iv]) / dy
    else
        return (w[ix, iy + 1, iv] - w[ix, iy - 1, iv]) / (2dy)
    end
end

@inline function _diff3_x(w, ix, iy, iz, nx, dx, iv)
    if ix == 1
        return (w[2, iy, iz, iv] - w[1, iy, iz, iv]) / dx
    elseif ix == nx
        return (w[nx, iy, iz, iv] - w[nx - 1, iy, iz, iv]) / dx
    else
        return (w[ix + 1, iy, iz, iv] - w[ix - 1, iy, iz, iv]) / (2dx)
    end
end

@inline function _diff3_y(w, ix, iy, iz, ny, dy, iv)
    if iy == 1
        return (w[ix, 2, iz, iv] - w[ix, 1, iz, iv]) / dy
    elseif iy == ny
        return (w[ix, ny, iz, iv] - w[ix, ny - 1, iz, iv]) / dy
    else
        return (w[ix, iy + 1, iz, iv] - w[ix, iy - 1, iz, iv]) / (2dy)
    end
end

@inline function _diff3_z(w, ix, iy, iz, nz, dz, iv)
    if iz == 1
        return (w[ix, iy, 2, iv] - w[ix, iy, 1, iv]) / dz
    elseif iz == nz
        return (w[ix, iy, nz, iv] - w[ix, iy, nz - 1, iv]) / dz
    else
        return (w[ix, iy, iz + 1, iv] - w[ix, iy, iz - 1, iv]) / (2dz)
    end
end
