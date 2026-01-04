using Batsrus
using PyPlot
using Printf
using Statistics
using GaussianMixtures
using Random
using LinearAlgebra
using StaticArrays
using FHist

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
   # We generate in local (para, perp1, perp2) and rotate
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
         # 80% Core, 20% Beam
         is_core = rand() > 0.2

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
            # We just need a vector perpendicular to B
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

# ---------------------
# Section 2: Coordinate Transformation
# ---------------------
println("\n--- Section 2: Coordinate Transformation ---")
B_field = [1.0, -1.0, 0.0]
println("Using B field: $B_field")

# Define dummy E field perpcendicular to B for full 3D basis
E_dummy = [1.0, 1.0, 0.0]
# Check orthogonality
if dot(B_field, E_dummy) != 0
   error("Dummy E field not perpendicular")
end

println("Using Dummy E field (for basis): $E_dummy")

# Get transformation function (returns v_para, v_perp1, v_perp2, weight)
transform_func = get_particle_field_aligned_transform(B_field, E_dummy)

# Apply transform manually for analysis
transformed_data, transformed_names = transform_func(
   particles, data.header.real_component_names)

v_para = transformed_data[1, :]
v_perp1 = transformed_data[2, :]
v_perp2 = transformed_data[3, :]
weights = transformed_data[4, :]

println("Transformed components: $transformed_names")

# ---------------------
# Section 3: Analysis & Plotting
# ---------------------
println("\n--- Section 3: Analysis & Plotting ---")

# Estimate Core parameters from histograms
# Helper to detect velocity names
function detect_velocity_names(rnames)
   vx_name = "vx"
   for cand in ["ux", "vx", "velocity_x"]
      if cand in rnames
         vx_name = cand
         break
      end
   end

   vy_name = "vy"
   for cand in ["uy", "vy", "velocity_y"]
      if cand in rnames
         vy_name = cand
         break
      end
   end
   return vx_name, vy_name
end

# Function to get peak and width from 1D weighted histogram
function estimate_1d_param(v, w, nbins = 50)
   h = Hist1D(v; binedges = range(minimum(v), maximum(v), length = nbins + 1), weights = w)
   _, max_idx = findmax(h.bincounts)
   edges = h.binedges isa Tuple ? h.binedges[1] : h.binedges
   peak_loc = (edges[max_idx] + edges[max_idx + 1]) / 2

   # Estimate width
   return peak_loc, 1.0
end

bulk_para, vth_para_est = estimate_1d_param(v_para, weights)
bulk_perp1, vth_perp1_est = estimate_1d_param(v_perp1, weights)
bulk_perp2, vth_perp2_est = estimate_1d_param(v_perp2, weights)

println("Estimated Core Bulk: Para=$bulk_para, Perp1=$bulk_perp1, Perp2=$bulk_perp2")

# Construct parameter vectors
bulk_vel = [bulk_para, bulk_perp1, bulk_perp2]
vth_vec = [vth_para_est, vth_perp1_est, vth_perp2_est]
vth_scalar = maximum(vth_vec)

# Classify
velocities_3d = transformed_data[1:3, :]

println("Classifying particles...")
mask_core = Batsrus.get_core_population_mask(velocities_3d, bulk_vel, vth_scalar, 3.0)

println("Core particles: ", count(mask_core))
println("Suprathermal particles: ", count(.!mask_core))

# Split data
core_data = transformed_data[:, mask_core]
supra_data = transformed_data[:, .!mask_core]

# Fit GMMs
# Function to fit and print
function fit_and_report(label, data_subset, n_comp)
   println("\nFitting $label (N=$(size(data_subset,2)))...")
   if size(data_subset, 2) < n_comp * 10
      println("  Not enough particles for GMM.")
      return nothing
   end

   # Extract velocities (3xN) and weights
   vels = data_subset[1:3, :]
   ws = data_subset[4, :]

   res = fit_particle_velocity_gmm(vels, n_comp; weights = ws)

   for (i, c) in enumerate(res)
      println("  Comp $i: W=$(c.weight), Mean=$(c.mean), Vth=$(c.vth)")
   end
   return res
end

res_core = fit_and_report("Core", core_data, 1)
res_supra = fit_and_report("Suprathermal", supra_data, 1) # Assuming single beam for now

# Combine results for reconstruction
all_results = []
total_weight = sum(weights)
weight_core = sum(core_data[4, :])
weight_supra = sum(supra_data[4, :])

if !isnothing(res_core)
   for c in res_core
      # Re-weight global weight
      push!(all_results,
         (
            weight = c.weight * (weight_core / total_weight),
            mean = c.mean,
            vth = c.vth
         ))
   end
end
if !isnothing(res_supra)
   for c in res_supra
      push!(all_results,
         (
            weight = c.weight * (weight_supra / total_weight),
            mean = c.mean,
            vth = c.vth
         ))
   end
end

fig, axs = plt.subplots(2, 2, figsize = (12, 12), constrained_layout = true)

# --- Plot 1: Original Coordinates (vx, vy) transformed via plot_phase ---

# Detect names
# Detect names
vx_name, vy_name = detect_velocity_names(data.header.real_component_names)

plot_phase(data, vx_name, vy_name, ax = axs[1],
   x_range = x_range, z_range = z_range,
   y_range = selection_y_range, bins = 100)
axs[1].set_title("Original Coordinates ($vx_name, $vy_name)")

# --- Plot 2: Transformed Coordinates (v_para, v_perp_mag) ---
# We need to construct v_perp magnitude from our 3D basis for the 2D plot.
v_perp_mag = sqrt.(v_perp1 .^ 2 .+ v_perp2 .^ 2)

h2 = Hist2D((v_para, v_perp_mag), nbins = (100, 100), weights = weights)
# Determine generic vmin/vmax from this main plot to use elsewhere
vmin_g = minimum(h2.bincounts[h2.bincounts .> 0])
vmax_g = maximum(h2.bincounts)

im2 = axs[2].imshow(h2.bincounts', origin = "lower",
   extent = [minimum(v_para), maximum(v_para), minimum(v_perp_mag), maximum(v_perp_mag)],
   norm = PyPlot.matplotlib.colors.LogNorm(vmin = vmin_g, vmax = vmax_g), aspect = "auto")
axs[2].set_xlabel("\$v_{\\parallel}\$")
axs[2].set_ylabel("\$v_{\\perp}\$")
axs[2].set_title("Field-Aligned Coordinates")
plt.colorbar(im2, ax = axs[2], pad = 0.02)

# --- Plot 3: GMM Reconstruction ---
axs[3].set_title("GMM Reconstruction (Split Fit)")
# Evaluate PDF on grid
xi = range(minimum(v_para), maximum(v_para), length = 100)
yi = range(minimum(v_perp_mag), maximum(v_perp_mag), length = 100)
Z = zeros(100, 100)

for i in 1:100, j in 1:100
   vp = xi[i]
   vt = yi[j] # v_perp magnitude

   prob = 0.0
   for res in all_results
      # res has 3D mean/vth: (para, perp1, perp2)
      mu_p = res.mean[1]
      sig_p = res.vth[1] / sqrt(2)

      # Average perp sigma for simplified Cylindrical PDF
      sig_t = sqrt((res.vth[2]^2 + res.vth[3]^2) / 4)

      # Parallel part
      p_para = exp(-0.5 * ((vp - mu_p) / sig_p)^2) / (sig_p * sqrt(2π))

      # Perp part (Rayleigh-like for magnitude)
      p_perp = (vt / sig_t^2) * exp(-0.5 * (vt / sig_t)^2)

      prob += res.weight * p_para * p_perp
   end
   # Scale prob to match histogram counts roughly (N * bin_area)
   # Or just plot PDF contours.
   Z[j, i] = prob
end

# Use histogram counts scale for comparison? No, just PDF contours is fine.
# But to match colorbar, we need to multiply by TotalWeight * BinArea.
# BinArea approx:
dx = step(xi)
dy = step(yi)
Z_scaled = Z .* (total_weight * dx * dy)

# Plot contours or heatmap
im3 = axs[3].imshow(Z_scaled, origin = "lower", extent = [xi[1], xi[end], yi[1], yi[end]],
   norm = PyPlot.matplotlib.colors.LogNorm(vmin = vmin_g, vmax = vmax_g), aspect = "auto")
plt.colorbar(im3, ax = axs[3], pad = 0.02)
axs[3].set_xlabel(L"v_{\parallel}")
axs[3].set_ylabel(L"v_{\perp}")

# --- Plot 4: Classification ---
axs[4].set_title("Classification")
# Scatter plot
# Subsample for speed if needed
stride = 1
if size(v_para, 1) > 10000
   stride = 10
end

axs[4].scatter(v_para[mask_core][1:stride:end], v_perp_mag[mask_core][1:stride:end],
   s = 1, c = "red", label = "Core", alpha = 0.5)
axs[4].scatter(v_para[.!mask_core][1:stride:end], v_perp_mag[.!mask_core][1:stride:end],
   s = 1, c = "blue", label = "Supra", alpha = 0.5)
axs[4].legend(markerscale = 5.0)
axs[4].set_xlabel(L"v_{\parallel}")
axs[4].set_ylabel(L"v_{\perp}")

savefig("phase_space_analysis.png")
println("Unified analysis plot saved to phase_space_analysis.png")
