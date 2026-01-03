using Batsrus
using PyPlot
using Printf
using Statistics
using GaussianMixtures
using Random
using LinearAlgebra
using StaticArrays

# Script to load AMReX particle data and plot phase space distribution
# Usage: julia examples/plot_phase_distribution.jl

# --- Configuration ---
# Note: Update this path to your specific environment
dir_prefix = "/global/homes/h/hyzhou/scratch/swmf/pleiades/2dmagnetosphere/run_2dXZ_f1_950s/PC/"
data_path = dir_prefix * "cut_particle_region0_5_t00001550_n00169931_amrex"

# Region of interest
x_center = 18.0
z_center = -5.0
width = 1.0

# ---------------------
# Section 1: Data Loading & Selection
# ---------------------
println("--- Section 1: Data Loading & Selection ---")

# Check if real data exists, otherwise generate mock data
if !isdir(data_path)
   println("Real data not found at: $data_path")
   println("Generating mock AMReX data for demonstration...")

   data_path = "data_mock_demo"
   mkpath(data_path)

   # Define Field
   B_field = [1.0, -1.0, 0.0]
   B_hat = normalize(B_field)

   # Rotation matrix to align Z with B for generating anisotropic beam
   # We want local_parallel -> B_hat
   # We can use Batsrus.getRotationMatrix if available or construct manually
   # Let's use a simple construction:
   # Z' = B_hat
   # Y' = Z x B_hat / |...| (if B not along Z)
   # X' = Y' x Z'

   # Or simpler: generate in local (para, perp1, perp2) and rotate
   function random_perp_vector(n)
      v = randn(3)
      v -= dot(v, n) * n
      return normalize(v)
   end

   Batsrus.generate_mock_amrex_data(
      data_path,
      num_particles = 10000,
      # Names according to typical AMReX output
      real_component_names = ["ux", "uy", "uz", "weight"],
      particle_gen = (i, n_reals) -> begin
         # 50% Core, 50% Beam
         is_core = rand() > 0.5

         if is_core
            # Core: Centered at (-1, 1, 0), Isotropic vth=1
            center = [-1.0, 1.0, 0.0]
            v_vec = center .+ randn(3)
         else
            # Beam: Centered at (2, -2, 0)
            # Anisotropic: vth_para = 0.5, vth_perp = 2.0
            center = [2.0, -2.0, 0.0]

            # Generate in local frame where z is parallel to B
            v_local_para = randn() * 0.5
            v_local_perp1 = randn() * 2.0
            v_local_perp2 = randn() * 2.0

            # Construct global velocity from local components
            # v_global = v_para * b_hat + v_perp1 * e1 + v_perp2 * e2
            # Since we don't care about specific perp direction orientation, 
            # we just need a vector perpendicular to B
            perp_vec = random_perp_vector(B_hat)
            perp_vec2 = cross(B_hat, perp_vec)

            v_vec = center .+ (v_local_para .* B_hat) .+ (v_local_perp1 .* perp_vec) .+
                    (v_local_perp2 .* perp_vec2)
         end

         return (
            0.5 + (rand() - 0.5) * 0.4, 0.5 + (rand() - 0.5) * 0.4, 0.5 +
                                                                    (rand() - 0.5) * 0.4, # x, y, z
            v_vec[1], v_vec[2], v_vec[3], # ux, uy, uz
            1.0 # weight
         )
      end
   )
   println("Mock data generated at: $data_path")

   # Update selection center for mock data
   global x_center = 0.5
   global z_center = 0.5
   global width = 0.6
end

# Load the data
data = AMReXParticle(data_path)

# Define ranges
half_width = width / 2
x_range = (x_center - half_width, x_center + half_width)
z_range = (z_center - half_width, z_center + half_width)
selection_y_range = (data.dim == 2) ? z_range : nothing

# Select particles
println("Selecting particles...")
if data.dim == 2
   particles = select_particles_in_region(data; x_range, y_range = selection_y_range)
else
   particles = select_particles_in_region(data; x_range, z_range)
end
println("Number of particles: ", size(particles, 2))

# Extract original velocities
names = data.header.real_component_names
idx_vx = findfirst(n -> n in ["ux", "vx", "velocity_x"], names)
idx_vy = findfirst(n -> n in ["uy", "vy", "velocity_y"], names)
idx_vz = findfirst(n -> n in ["uz", "vz", "velocity_z"], names)

vx = particles[idx_vx, :]
vy = particles[idx_vy, :]
vz = (idx_vz !== nothing) ? particles[idx_vz, :] : zeros(length(vx))

# ---------------------
# Section 2: Coordinate Transformation
# ---------------------
println("\n--- Section 2: Coordinate Transformation ---")
B_field = [1.0, -1.0, 0.0]
println("Using B field: $B_field")

# Get transformation function
transform_func = get_particle_field_aligned_transform(B_field)

# Apply transform manually for plotting/analysis matrix
# Typically this returns (new_data, new_names)
# passing select_particles_in_region result (which contains ALL real components)
transformed_data, transformed_names = transform_func(
   particles, data.header.real_component_names)

v_para = transformed_data[1, :]
v_perp = transformed_data[2, :]

println("Transformed components: $transformed_names")

# ---------------------
# Section 3: Plotting & GMM
# ---------------------
println("\n--- Section 3: Plotting & GMM ---")

fig, axs = plt.subplots(2, 2, figsize = (12, 12), constrained_layout = true)

# --- Plot 1: Original Coordinates (vx, vy) ---
ax1 = axs[1]
h1 = ax1.hist2d(
   vx, vy, bins = 100, norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1), cmap = "viridis")
ax1.set_xlabel("\$v_x\$")
ax1.set_ylabel("\$v_y\$")
ax1.set_title("Original Coordinates (vx, vy)")
plt.colorbar(h1[4], ax = ax1)

# --- Plot 2: Transformed Coordinates (v_para, v_perp) ---
ax2 = axs[2]
h2 = ax2.hist2d(v_para, v_perp, bins = 100,
   norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1), cmap = "viridis")
ax2.set_xlabel("\$v_{\\parallel}\$")
ax2.set_ylabel("\$v_{\\perp}\$")
ax2.set_title("Field-Aligned Coordinates")
plt.colorbar(h2[4], ax = ax2)

# --- GMM Fitting on Transformed Data ---
println("Fitting GMM on transformed data...")
n_clusters = 2
# Pass transformed matrix (2 x N) to the new matrix-input dispatch
# Note: transformed_data is (2, N)
gmm_results = fit_particle_velocity_gmm(transformed_data, n_clusters)

println("GMM Results (Transformed Frame):")
for (i, res) in enumerate(gmm_results)
   println("  Component $i: Weight=$(res.weight), Mean=$(res.mean), Vth=$(res.vth)")
end

# --- Plot 3: GMM Reconstruction (Transformed) ---
ax3 = axs[3]
# Plot histogram again as background
ax3.hist2d(v_para, v_perp, bins = 100,
   norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1), cmap = "Greys")

# Overlay GMM contours
# Grid for evaluation
x_min, x_max = minimum(v_para), maximum(v_para)
y_min, y_max = minimum(v_perp), maximum(v_perp)
xi = range(x_min, x_max, length = 100)
yi = range(y_min, y_max, length = 100)
Z = zeros(100, 100)

for i in 1:100, j in 1:100
   v_loc = [xi[i], yi[j]] # (v_para, v_perp)
   prob = 0.0
   for k in 1:n_clusters
      # Diagonal covariance assumption in fit_particle_velocity_gmm (kind=:diag)
      # P(v) = w * PDF(v_para) * PDF(v_perp)
      mu = gmm_results[k].mean
      sigma = gmm_results[k].vth ./ sqrt(2)

      p_para = exp(-0.5 * ((v_loc[1] - mu[1]) / sigma[1])^2) / (sigma[1] * sqrt(2π))
      p_perp = exp(-0.5 * ((v_loc[2] - mu[2]) / sigma[2])^2) / (sigma[2] * sqrt(2π))

      prob += gmm_results[k].weight * p_para * p_perp
   end
   Z[j, i] = prob
end

levels = exp.(range(log(maximum(Z) * 1e-4), log(maximum(Z)), length = 10))
ax3.contour(xi, yi, Z, levels = levels, cmap = "viridis", linewidths = 2)
ax3.set_xlabel("\$v_{\\parallel}\$")
ax3.set_ylabel("\$v_{\\perp}\$")
ax3.set_title("GMM Reconstruction (Transformed)")

# --- Plot 4: Classification (Transformed) ---
ax4 = axs[4]
# Assign particles to clusters (Hard assignment for visualization)
labels = Vector{Int}(undef, length(v_para))
for i in eachindex(v_para)
   best_k = 0
   max_p = -1.0
   v_vec = [v_para[i], v_perp[i]]

   for k in 1:n_clusters
      mu = gmm_results[k].mean
      sigma = gmm_results[k].vth ./ sqrt(2)
      p = gmm_results[k].weight *
          (exp(-0.5 * ((v_vec[1] - mu[1]) / sigma[1])^2) / sigma[1]) *
          (exp(-0.5 * ((v_vec[2] - mu[2]) / sigma[2])^2) / sigma[2])
      if p > max_p
         max_p = p
         best_k = k
      end
   end
   labels[i] = best_k
end

# Color by label
colors = ["red", "blue"]
for k in 1:n_clusters
   mask = labels .== k
   if any(mask)
      ax4.scatter(v_para[mask], v_perp[mask], s = 1, c = colors[k],
         label = "Cluster $k", alpha = 0.5)
   end
end
ax4.set_xlabel("\$v_{\\parallel}\$")
ax4.set_ylabel("\$v_{\\perp}\$")
ax4.set_title("GMM Classification")
ax4.legend()

savefig("phase_space_analysis.png")
println("Unified analysis plot saved to phase_space_analysis.png")
