# Derived quantities from raw output variables.

# Earth radius in km used for current density scaling in PLANETARY units.
const EARTH_RADIUS_KM = 6378.0
# Elementary charge in C
const ELEMENTARY_CHARGE = 1.602176634e-19
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
    return _get_magnitude(bd, Val(var))
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
            res[i] = sqrt(
                w[i + (ivx - 1) * n_space]^2 +
                    w[i + (ivy - 1) * n_space]^2 +
                    w[i + (ivz - 1) * n_space]^2
            )
        end
        return rebuild(bd.w, res, (d[1], d[2]))
    elseif ndim == 3
        nx, ny, nz = size(bd.w, 1), size(bd.w, 2), size(bd.w, 3)
        n_space = nx * ny * nz
        res = similar(w, nx, ny, nz)
        @inbounds for i in 1:n_space
            res[i] = sqrt(
                w[i + (ivx - 1) * n_space]^2 +
                    w[i + (ivy - 1) * n_space]^2 +
                    w[i + (ivz - 1) * n_space]^2
            )
        end
        return rebuild(bd.w, res, (d[1], d[2], d[3]))
    else
        sz = ntuple(i -> size(bd.w, i), Val(ndim))
        n_space = 1
        for i in 1:ndim
            n_space *= sz[i]
        end
        res = similar(w, sz)
        @inbounds for i in 1:n_space
            res[i] = sqrt(
                w[i + (ivx - 1) * n_space]^2 +
                    w[i + (ivy - 1) * n_space]^2 +
                    w[i + (ivz - 1) * n_space]^2
            )
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
        res[i] = sqrt(
            w[i + (ivx - 1) * n_cells]^2 +
                w[i + (ivy - 1) * n_cells]^2 +
                w[i + (ivz - 1) * n_cells]^2
        )
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

@inline function _get_density_index(bd::BatsrusIDL, species)
    return findindex(bd, "rhos" * string(species))
end

@inline function _get_species_mass(bd::BatsrusIDL, species)
    # Check for species mass in parameters
    for prefix in ("ms", "mass")
        idx = findlast(x -> x == prefix * string(species), bd.head.param)
        if !isnothing(idx)
            return bd.head.eqpar[idx]
        end
    end
    # Fallback to 1.0 and warn if it's electrons (species 0)
    if species == 0
        @warn "Species 0 mass not found in parameters; using 1.0"
    end
    return 1.0
end

@inline function _get_pe_scaling(bd::BatsrusIDL)
    if startswith(bd.head.headline, "normalized")
        return 1.0
    elseif bd.head.headline == "PLANETARY" || occursin(" nT ", bd.head.headline)
        # Factor to convert (nPa / Re) / (amu/cc * e) to μV/m
        return 1.0e-12 / (EARTH_RADIUS_KM * ELEMENTARY_CHARGE)
    else
        # Factor to convert (nPa / km) / (amu/cc * e) to μV/m
        return 1.0e-12 / ELEMENTARY_CHARGE
    end
end

@inline function _get_hall_scaling(bd::BatsrusIDL, hasJ::Bool)
    if startswith(bd.head.headline, "normalized")
        return 1.0
    elseif hasJ || bd.head.headline == "PLANETARY"
        # Factor to convert (μA/m² * nT) / (amu/cc) to μV/m
        # E [V/m] = (10⁻⁶ * 10⁻⁹) / (10⁶ * e) * (j*b/n) = 10⁻²¹ / e * (j*b/n)
        # E [μV/m] = 10⁻¹⁵ / e * (j*b/n)
        return 1.0e-15 / ELEMENTARY_CHARGE
    else
        # Factor to convert (raw_J * nT) / (amu/cc) to μV/m
        # J_SI = (raw_J * 10⁻⁹ / 10³) / μ₀ = raw_J * 10⁻¹² / (4π*10⁻⁷)
        # E [V/m] = (raw_J * 10⁻¹² / μ₀ * 10⁻⁹) / (10⁶ * e) = 10⁻²⁷ / (μ₀*e) * (j*b/n)
        # E [μV/m] = 10⁻²¹ / (μ₀*e) * (j*b/n)
        return 1.0e-21 / (4.0 * π * 1.0e-7 * ELEMENTARY_CHARGE)
    end
end


"""
    get_anisotropy(bd::BatsrusIDL, species=0; method=:projection)

Calculate the pressure anisotropy for `species`. Two methods are supported:
- `:projection`: direct projection of the pressure tensor onto the magnetic field (default).
- `:rotation`: rotating the pressure tensor to a field-aligned coordinate system.
"""
@inline function get_anisotropy(bd::BatsrusIDL{ndim, TV}, species = 0; method = :projection) where {ndim, TV}
    return _get_anisotropy(bd, species, method)
end

@inline function _get_B_P(w, i, n, iv, ip)
    ivx, ivy, ivz = iv
    ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip
    Bx = w[i + (ivx - 1) * n]
    By = w[i + (ivy - 1) * n]
    Bz = w[i + (ivz - 1) * n]
    pxx = w[i + (ipxx - 1) * n]
    pyy = w[i + (ipyy - 1) * n]
    pzz = w[i + (ipzz - 1) * n]
    pxy = w[i + (ipxy - 1) * n]
    pxz = w[i + (ipxz - 1) * n]
    pyz = w[i + (ipyz - 1) * n]
    return Bx, By, Bz, pxx, pyy, pzz, pxy, pxz, pyz
end

@inline function _calc_anisotropy(Bx, By, Bz, pxx, pyy, pzz, pxy, pxz, pyz, method)
    if method === :projection
        B2 = Bx^2 + By^2 + Bz^2
        p_parallel = (pxx * Bx^2 + pyy * By^2 + pzz * Bz^2 + 2 * pxy * Bx * By + 2 * pxz * Bx * Bz + 2 * pyz * By * Bz) / B2
        return (pxx + pyy + pzz - p_parallel) / (2 * p_parallel)
    elseif method === :rotation
        P = @SMatrix [pxx pxy pxz; pxy pyy pyz; pxz pyz pzz]
        v = @SVector [Bx, By, Bz]
        Pr = rotateTensorToVectorZ(P, v)
        return (Pr[1, 1] + Pr[2, 2]) / (2 * Pr[3, 3])
    else
        error("Unknown method $method")
    end
end

@inline function _get_anisotropy(bd::BatsrusIDLStructured{ndim, TV, TX, TW}, species, method) where {ndim, TV, TX, TW}
    iv = get_vectors_indices(bd, :B)
    ip = _get_pressure_tensor_indices(bd, species)

    w = parent(bd.w)
    d = dims(bd.w)

    sz = ntuple(i -> size(bd.w, i), Val(ndim))
    n_space = prod(sz)
    res = similar(w, sz)

    @inbounds for i in 1:n_space
        vars = _get_B_P(w, i, n_space, iv, ip)
        res[i] = _calc_anisotropy(vars..., method)
    end
    return rebuild(bd.w, res, ntuple(i -> d[i], Val(ndim)))
end

@inline function _get_anisotropy(bd::BatsrusIDLUnstructured{ndim, TV, TX, TW}, species, method) where {ndim, TV, TX, TW}
    iv = get_vectors_indices(bd, :B)
    ip = _get_pressure_tensor_indices(bd, species)

    w = parent(bd.w)
    n_cells = size(w, 1)

    res = similar(w, n_cells)

    @inbounds for i in 1:n_cells
        vars = _get_B_P(w, i, n_cells, iv, ip)
        res[i] = _calc_anisotropy(vars..., method)
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
function get_hall_E(bd::BatsrusIDL{ndim, TV}) where {ndim, TV}
    Bx, By, Bz = get_vectors(bd, :B)

    if _has_var(bd, "uxs0")
        uex, uey, uez = get_vectors(bd, :U0)
        # Let us use H+ velocities as the ion bulk velocity and ignore heavy ions
        uix, uiy, uiz = get_vectors(bd, :U1)

        Ehallx = similar(Bx)
        Ehally = similar(By)
        Ehallz = similar(Bz)
        # (Ui - Ue) × B
        @inbounds @simd for i in eachindex(Ehallx)
            Ehallx[i] = (uiy[i] - uey[i]) * Bz[i] - (uiz[i] - uez[i]) * By[i]
            Ehally[i] = (uiz[i] - uez[i]) * Bx[i] - (uix[i] - uex[i]) * Bz[i]
            Ehallz[i] = (uix[i] - uex[i]) * By[i] - (uiy[i] - uey[i]) * Bx[i]
        end
    else
        hasJ = _has_var(bd, "jx")
        Jx, Jy, Jz = get_current_density(bd)
        rho = _has_var(bd, "rho") ? getvar(bd, "rho") : getvar(bd, "rhos1")
        # Ions are species 1 by convention
        m = _get_species_mass(bd, 1)
        C = TV(_get_hall_scaling(bd, hasJ))

        Ehallx = similar(Bx)
        Ehally = similar(By)
        Ehallz = similar(Bz)

        # E = C * J × B / (rho/m)
        @inbounds @simd for i in eachindex(Ehallx)
            ne_inv = m / rho[i]
            Ehallx[i] = C * ne_inv * (Jy[i] * Bz[i] - Jz[i] * By[i])
            Ehally[i] = C * ne_inv * (Jz[i] * Bx[i] - Jx[i] * Bz[i])
            Ehallz[i] = C * ne_inv * (Jx[i] * By[i] - Jy[i] * Bx[i])
        end
    end

    return Ehallx, Ehally, Ehallz
end

"""
    get_pe_E(bd::BatsrusIDL, species=0; mass=nothing)

Return the electric field ``\\mathbf{E}_{p_e} = -\\frac{1}{n_e e} \\nabla \\cdot \\mathbf{P}_e``
derived from the divergence of the electron pressure tensor. Units are in μV/m if
PLANETARY or km units are used, otherwise 1.0.
Note that species 0 is by convention electrons, but it is not guaranteed.
"""
@inline function get_pe_E(bd::BatsrusIDL, species = 0; mass = nothing)
    return _get_pe_E(bd, species, mass)
end

function _get_pe_E(bd::BatsrusIDLStructured{1, TV}, species, mass) where {TV}
    irho = _get_density_index(bd, species)
    ip = _get_pressure_tensor_indices(bd, species)
    ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip

    m = isnothing(mass) ? _get_species_mass(bd, species) : TV(mass)
    C = TV(_get_pe_scaling(bd))

    w = parent(bd.w)
    xrange = get_range(bd)[1]
    dx, nx = TV(step(xrange)), size(w, 1)

    ex = zeros(TV, nx)
    ey = zeros(TV, nx)
    ez = zeros(TV, nx)

    for i in 1:nx
        rho = w[i, irho]
        ne = rho / m
        # Div P = [d_x Pxx, d_x Pxy, d_x Pxz]
        ex[i] = -C / ne * _diff1(w, i, nx, dx, ipxx)
        ey[i] = -C / ne * _diff1(w, i, nx, dx, ipxy)
        ez[i] = -C / ne * _diff1(w, i, nx, dx, ipxz)
    end

    dims_ = dims(bd.w)[1:1]
    return DimArray(ex, dims_), DimArray(ey, dims_), DimArray(ez, dims_)
end

function _get_pe_E(bd::BatsrusIDLStructured{2, TV}, species, mass) where {TV}
    irho = _get_density_index(bd, species)
    ip = _get_pressure_tensor_indices(bd, species)
    ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip

    m = isnothing(mass) ? _get_species_mass(bd, species) : TV(mass)
    C = TV(_get_pe_scaling(bd))

    w = parent(bd.w)
    xrange, yrange = get_range(bd)
    dx, dy = TV(step(xrange)), TV(step(yrange))
    nx, ny = size(w, 1), size(w, 2)

    ex = zeros(TV, nx, ny)
    ey = zeros(TV, nx, ny)
    ez = zeros(TV, nx, ny)

    for iy in 1:ny, ix in 1:nx
        rho = w[ix, iy, irho]
        ne = rho / m
        # (Div P)_x = d_x Pxx + d_y Pxy
        dpxx_dx = _diff2_x(w, ix, iy, nx, dx, ipxx)
        dpxy_dy = _diff2_y(w, ix, iy, ny, dy, ipxy)
        ex[ix, iy] = -C / ne * (dpxx_dx + dpxy_dy)

        # (Div P)_y = d_x Pxy + d_y Pyy
        dpxy_dx = _diff2_x(w, ix, iy, nx, dx, ipxy)
        dpyy_dy = _diff2_y(w, ix, iy, ny, dy, ipyy)
        ey[ix, iy] = -C / ne * (dpxy_dx + dpyy_dy)

        # (Div P)_z = d_x Pxz + d_y Pyz
        dpxz_dx = _diff2_x(w, ix, iy, nx, dx, ipxz)
        dpyz_dy = _diff2_y(w, ix, iy, ny, dy, ipyz)
        ez[ix, iy] = -C / ne * (dpxz_dx + dpyz_dy)
    end

    dims_ = dims(bd.w)[1:2]
    return DimArray(ex, dims_), DimArray(ey, dims_), DimArray(ez, dims_)
end

function _get_pe_E(bd::BatsrusIDLStructured{3, TV}, species, mass) where {TV}
    irho = _get_density_index(bd, species)
    ip = _get_pressure_tensor_indices(bd, species)
    ipxx, ipyy, ipzz, ipxy, ipxz, ipyz = ip

    m = isnothing(mass) ? _get_species_mass(bd, species) : TV(mass)
    C = TV(_get_pe_scaling(bd))

    w = parent(bd.w)
    xrange, yrange, zrange = get_range(bd)
    dx, dy, dz = TV(step(xrange)), TV(step(yrange)), TV(step(zrange))
    nx, ny, nz = size(w, 1), size(w, 2), size(w, 3)

    ex = zeros(TV, nx, ny, nz)
    ey = zeros(TV, nx, ny, nz)
    ez = zeros(TV, nx, ny, nz)

    for iz in 1:nz, iy in 1:ny, ix in 1:nx
        rho = w[ix, iy, iz, irho]
        ne = rho / m

        # (Div P)_x = d_x Pxx + d_y Pxy + d_z Pxz
        dpxx_dx = _diff3_x(w, ix, iy, iz, nx, dx, ipxx)
        dpxy_dy = _diff3_y(w, ix, iy, iz, ny, dy, ipxy)
        dpxz_dz = _diff3_z(w, ix, iy, iz, nz, dz, ipxz)
        ex[ix, iy, iz] = -C / ne * (dpxx_dx + dpxy_dy + dpxz_dz)

        # (Div P)_y = d_x Pxy + d_y Pyy + d_z Pyz
        dpxy_dx = _diff3_x(w, ix, iy, iz, nx, dx, ipxy)
        dpyy_dy = _diff3_y(w, ix, iy, iz, ny, dy, ipyy)
        dpyz_dz = _diff3_z(w, ix, iy, iz, nz, dz, ipyz)
        ey[ix, iy, iz] = -C / ne * (dpxy_dx + dpyy_dy + dpyz_dz)

        # (Div P)_z = d_x Pxz + d_y Pyz + d_z Pzz
        dpxz_dx = _diff3_x(w, ix, iy, iz, nx, dx, ipxz)
        dpyz_dy = _diff3_y(w, ix, iy, iz, ny, dy, ipyz)
        dpzz_dz = _diff3_z(w, ix, iy, iz, nz, dz, ipzz)
        ez[ix, iy, iz] = -C / ne * (dpxz_dx + dpyz_dy + dpzz_dz)
    end

    dims_ = dims(bd.w)[1:3]
    return DimArray(ex, dims_), DimArray(ey, dims_), DimArray(ez, dims_)
end

function _get_pe_E(bd::BatsrusIDLUnstructured, species, mass)
    error("get_pe_E is not yet supported for unstructured grids.")
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
    _has_var(bd, "jx") && return getvar(bd, "jx")
    return _compute_jx(bd)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:jy})
    _has_var(bd, "jy") && return getvar(bd, "jy")
    return _compute_jy(bd)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:jz})
    _has_var(bd, "jz") && return getvar(bd, "jz")
    return _compute_jz(bd)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:j})
    if _has_var(bd, "j")
        return getvar(bd, "j")
    elseif _has_var(bd, "jx") && _has_var(bd, "jy") && _has_var(bd, "jz")
        # Reuse get_magnitude if indices can be found
        return get_magnitude(bd, :J)
    end
    Jx, Jy, Jz = get_current_density(bd)
    return sqrt.(Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2)
end

@inline function _getvar(bd::BatsrusIDL, ::Val{:b})
    _has_var(bd, "b") && return getvar(bd, "b")
    return get_magnitude(bd, :B)
end
@inline function _getvar(bd::BatsrusIDL, ::Val{:b2})
    _has_var(bd, "b2") && return getvar(bd, "b2")
    return get_magnitude2(bd, :B)
end
@inline function _getvar(bd::BatsrusIDL, ::Val{:e})
    _has_var(bd, "e") && return getvar(bd, "e")
    return get_magnitude(bd, :E)
end
@inline function _getvar(bd::BatsrusIDL, ::Val{:u})
    _has_var(bd, "u") && return getvar(bd, "u")
    return get_magnitude(bd, :U)
end
@inline _getvar(bd::BatsrusIDL, ::Val{:anisotropy0}) = get_anisotropy(bd, 0)
@inline _getvar(bd::BatsrusIDL, ::Val{:anisotropy1}) = get_anisotropy(bd, 1)

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
    jy = [-_diff1(w, ix, nx, dx, ivz) for ix in 1:nx]
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
    d = dims(bd.w)
    dx, dy = TV(step(val(d[1]))), TV(step(val(d[2])))
    nx, ny = size(w, 1), size(w, 2)
    jz = [
        _diff2_x(w, ix, iy, nx, dx, ivy) - _diff2_y(w, ix, iy, ny, dy, ivx)
            for ix in 1:nx, iy in 1:ny
    ]
    _apply_j_scaling!(jz, bd)
    return DimArray(jz, (d[1], d[2]))
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
