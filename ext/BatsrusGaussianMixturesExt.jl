module BatsrusGaussianMixturesExt

using Batsrus
import Batsrus: fit_particle_velocity_gmm, get_gmm_thermal_velocity
using GaussianMixtures: GMM, covars
using LinearAlgebra: diag, diagm, normalize, tr, dot
using StaticArrays: SVector
using PrecompileTools: @setup_workload, @compile_workload

"""
    fit_particle_velocity_gmm(data, n_clusters; x_range=nothing, y_range=nothing, z_range=nothing, vdim=3)

Fit a Gaussian Mixture Model to particle velocities in a region.

# Arguments

  - `data`: AMReXParticle data.
  - `n_clusters`: Number of GMM components.
  - `vdim`: Velocity dimension (1, 2, or 3).
  - `kind`: Covariance kind, `:full` (default) or `:diag`.

# Returns

  - A vector of named tuples sorted by weight, each containing:

      + `weight`: Component weight.
      + `mean`: Component mean velocity (vector of length vdim).
      + `cov`: Component covariance matrix (vdim x vdim).
      + `vth`: Component thermal velocity (diagonal approximation).
"""
function Batsrus.fit_particle_velocity_gmm(
        data::AMReXParticle{T},
        n_clusters::Int;
        x_range = nothing,
        y_range = nothing,
        z_range = nothing,
        vdim::Int = 3,
        kind::Symbol = :full
    ) where {T}
    particles = select_particles_in_region(data; x_range, y_range, z_range)
    if isempty(particles)
        error("No particles found in the specified region.")
    end

    vel_indices = Batsrus._get_velocity_indices(data, vdim)
    velocities = particles[vel_indices, :]

    return fit_particle_velocity_gmm(velocities, n_clusters; kind)
end

"""
    fit_particle_velocity_gmm(velocities::AbstractMatrix, n_clusters::Int; weights=nothing, kind=:full)

Fit a Gaussian Mixture Model to particle velocities.
"""
function Batsrus.fit_particle_velocity_gmm(
        velocities::AbstractMatrix{T},
        n_clusters::Int;
        weights::Union{AbstractVector, Nothing} = nothing,
        kind::Symbol = :full
    ) where {T}
    # GaussianMixtures expects (n_samples, n_features).
    vdim = size(velocities, 1)

    if isnothing(weights)
        X = Matrix{T}(velocities')
    else
        n_particles = size(velocities, 2)
        w_sum = sum(weights)
        probs = weights ./ w_sum
        cdf = cumsum(probs)
        cdf[end] = 1.0

        X = Matrix{T}(undef, n_particles, vdim)
        for i in 1:n_particles
            r = rand()
            idx = searchsortedfirst(cdf, r)
            idx = clamp(idx, 1, n_particles)
            X[i, :] = velocities[:, idx]
        end
    end

    gmm = GMM(n_clusters, X, kind = kind)

    results = map(1:n_clusters) do i
        μ = T.(gmm.μ[i, :])
        w = T(gmm.w[i])

        Σ = if kind == :diag
            diagm(T.(gmm.Σ[i, :]))
        else # :full
            T.(covars(gmm)[i])
        end

        vth_diag = sqrt.(2 .* diag(Σ))

        (weight = w, mean = μ, cov = Σ, vth = vth_diag)
    end

    sort!(results, by = x -> x.weight, rev = true)

    return results
end

"""
    get_gmm_thermal_velocity(gmm_component, b_field)

Calculate parallel and perpendicular thermal velocities from a GMM component.

# Returns

  - `(v_th_para, v_th_perp)`
"""
function Batsrus.get_gmm_thermal_velocity(
        gmm_component, b_field::AbstractVector{T}
    ) where {T}
    Σ = gmm_component.cov

    bhat = SVector{3, T}(normalize(b_field))

    var_para = dot(bhat, Σ * bhat)

    # Variance perpendicular to B (assuming gyrotropy)
    # Trace is invariant: tr(Σ) = var_para + 2 * var_perp
    var_perp = (tr(Σ) - var_para) / 2

    v_th_para = sqrt(2 * var_para)
    v_th_perp = sqrt(2 * var_perp)

    return v_th_para, v_th_perp
end

# Precompilation
@setup_workload begin
    @compile_workload begin
        mktempdir() do tmpdir
            Batsrus.generate_mock_amrex_data(
                tmpdir;
                real_component_names = ["ux", "uy", "uz"],
                particle_gen = (i, n_reals) -> (
                    rand(), rand(), rand(), randn(), randn(), randn(),
                )
            )
            data = AMReXParticle(tmpdir)

            fit_particle_velocity_gmm(data, 1; kind = :full)

            particles = select_particles_in_region(data)
            transform_func =
                Batsrus.get_particle_field_aligned_transform([1.0, 0.0, 0.0])
            transformed_data, _ = transform_func(
                particles, data.header.real_component_names
            )
            vels = transformed_data[1:3, :]

            gmm = fit_particle_velocity_gmm(vels, 1; kind = :diag)
            Batsrus.get_gmm_thermal_velocity(gmm[1], [1.0, 0.0, 0.0])
        end
    end
end

end
