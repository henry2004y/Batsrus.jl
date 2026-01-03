using Batsrus
using PyPlot
using Printf
using Statistics
using GaussianMixtures
using Random

# Script to load AMReX particle data and plot phase space distribution
# Usage: julia examples/plot_phase_distribution.jl

# --- Configuration ---
# Note: Update this path to your specific environment
dir_prefix = "/global/homes/h/hyzhou/scratch/swmf/pleiades/2dmagnetosphere/run_2dXZ_f1_950s/PC/"
data_path = dir_prefix * "cut_particle_region0_5_t00001550_n00169931_amrex"

# Region of interest
# Box centered at x=18, z=-5 with widths 1
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
   # Clean up previous mock data if exists to ensure fresh start or just reuse
   mkpath(data_path)

   # Generate mock data around 0.5 (assuming default domain is [0,1] or similar, 
   # or just ensure we match the selection range).
   # Let's generate in [0, 1] essentially.
   Batsrus.generate_mock_amrex_data(
      data_path,
      num_particles = 10000,
      real_component_names = ["ux", "uy", "uz", "weight"],
      particle_gen = (i, n_reals) -> begin
         # Create a mix of core and suprathermal particles with shifted centers
         is_core = rand() > 0.3 # 70% Core, 30% Beam/Halo
         if is_core
            # Core: Centered at 0, thermal velocity ~1
            vx = randn()
            vy = randn()
            vz = randn()
         else
            # Suprathermal/Beam: Centered at 5 in vx, thermal velocity ~3, anisotropic
            vx = 5.0 + randn() * 3.0
            vy = 3.0 * randn()
            vz = 3.0 * randn()
         end

         return (
            0.5 + (rand() - 0.5) * 0.4, 0.5 + (rand() - 0.5) * 0.4, 0.5 +
                                                                    (rand() - 0.5) * 0.4, # x, y, z near 0.5
            vx, vy, vz, # ux, uy, uz
            1.0 # weight
         )
      end
   )
   println("Mock data generated at: $data_path")

   # Update selection center for mock data
   global x_center = 0.5
   global z_center = 0.5
   global width = 0.6 # Cover the generated range
end

# Load the data
data = AMReXParticle(data_path)

# Define ranges
half_width = width / 2
x_range = (x_center - half_width, x_center + half_width)
z_range = (z_center - half_width, z_center + half_width)

println("Selecting particles in region:")
println("  x: $x_range")
println("  z: $z_range")

# Handle dimension mapping for 2D X-Z simulations
if data.dim == 2
   println("Detected 2D data. Assuming X-Z simulation, mapping z_range to y_range.")
   selection_y_range = z_range
else
   selection_y_range = nothing
end

# Select particles
if data.dim == 2
   particles = select_particles_in_region(
      data; x_range = x_range, y_range = selection_y_range)
else
   particles = select_particles_in_region(data; x_range = x_range, z_range = z_range)
end

println("Number of particles selected: ", size(particles, 2))

# Identify velocity components
names = data.header.real_component_names
function find_component_index(names, target_aliases)
   for (i, name) in enumerate(names)
      if name in target_aliases
         return i
      end
   end
   # Fallback for mock data standard names if not found in aliases
   for (i, name) in enumerate(names)
      if startswith(name, "u") && occursin(target_aliases[1][end:end], name)
         return i
      end
   end
   error("Component $(target_aliases) not found in data.")
end

idx_vx = find_component_index(names, ["vx", "velocity_x", "ux"])
idx_vy = find_component_index(names, ["vy", "velocity_y", "uy"])
# Try to find vz, might fail if 2D but mock data has it
idx_vz = try
   find_component_index(names, ["vz", "velocity_z", "uz"])
catch
   0
end

vx = particles[idx_vx, :]
vy = particles[idx_vy, :]
vz = idx_vz > 0 ? particles[idx_vz, :] : zeros(length(vx))

# ---------------------
# Unified Plotting & Analysis
# ---------------------
println("\n--- Unified Plotting & Analysis ---")

# Create a unified figure with 4 subplots
fig, axs = plt.subplots(2, 2, figsize = (14, 12), constrained_layout = true)

# ---------------------
# Subplot 1: Phase Space Distribution (Whole)
# ---------------------
println("Plotting Phase Space Distribution...")
ax1 = axs[1] # Top-Left
# Use plot_phase from Batsrus
plot_phase(
   data, names[idx_vx], names[idx_vy];
   x_range = x_range,
   y_range = selection_y_range, # Handle 2D mapping
   plot_zero_lines = true,
   edges = (range(-10, 10, length = 100), range(-10, 10, length = 100)),
   ax = ax1
)
ax1.set_title("Overall Distribution")

# ---------------------
# Subplot 2: GMM Fitting & Reconstruction
# ---------------------
println("Fitting and Plotting GMM...")
ax2 = axs[2] # Top-Right

println("Fitting Gaussian Mixture Model (k=2) via fit_particle_velocity_gmm...")

# Fit GMM using the new API
# Note: we use ranges to select data again (or could pass prepared data if API allowed)
results = fit_particle_velocity_gmm(
   data, 2;
   x_range = x_range,
   y_range = selection_y_range,
   vdim = 3
)

println("GMM Fit Results:")
for (i, comp) in enumerate(results)
   println("  Component $i:")
   println("    Weight: ", comp.weight)
   println("    Mean:   ", comp.mean)
   println("    Vth:    ", comp.vth)
end

# --- Plotting Reconstruction ---
# Plot original data histogram (Greys)
ax2.hist2d(
   vx, vy, bins = 100, norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1), cmap = "Greys"
)

# Create grid for contour plot
xmin, xmax = minimum(vx), maximum(vx)
ymin, ymax = minimum(vy), maximum(vy)

# Grid resolution
ngrid = 100
xi = range(xmin, xmax, length = ngrid)
yi = range(ymin, ymax, length = ngrid)

Z = zeros(ngrid, ngrid)

# Collect params for convenience
weights = [r.weight for r in results]
means = [r.mean for r in results]
vths = [r.vth for r in results] # vth = sqrt(2*sigma^2) => sigma = vth / sqrt(2)
sigmas = [v ./ sqrt(2) for v in vths]

gaussian_pdf(x, mu, sigma) = exp(-0.5 * ((x - mu) / sigma)^2) / (sigma * sqrt(2 * pi))

for i in 1:ngrid
   for j in 1:ngrid
      val = 0.0
      x_val = xi[i] # vx
      y_val = yi[j] # vy

      for k in 1:length(weights)
         # Marginal PDF(vx, vy)
         pdf_x = gaussian_pdf(x_val, means[k][1], sigmas[k][1])
         pdf_y = gaussian_pdf(y_val, means[k][2], sigmas[k][2])
         val += weights[k] * pdf_x * pdf_y
      end
      Z[j, i] = val # Note: PyPlot expects Z[row, col] -> y, x
   end
end

# Overlay contours
levels = exp.(range(log(maximum(Z) * 1e-4), log(maximum(Z)), length = 10))
c = ax2.contour(xi, yi, Z, levels = levels, cmap = "viridis", linewidths = 2)
ax2.clabel(c, inline = 1, fontsize = 8)

ax2.set_xlabel("\$v_x\$")
ax2.set_ylabel("\$v_y\$")
ax2.set_title("GMM Reconstruction")

# ---------------------
# Subplots 3 & 4: Classification
# ---------------------
println("Classifying Particles...")

# Estimate thermal velocity
vth_x = std(vx)
vth_y = std(vy)
vth_z = std(vz)
vth_est = sqrt((vth_x^2 + vth_y^2 + vth_z^2) / 3)
println("Estimated thermal velocity (vth): ", vth_est)

core, supra = classify_particles(
   data;
   x_range = x_range,
   y_range = selection_y_range,
   vdim = 3,
   vth = vth_est,
   nsigma = 3.0
)

println("Core particles: ", size(core, 2))
println("Suprathermal particles: ", size(supra, 2))

v_range = [[minimum(vx), maximum(vx)], [minimum(vy), maximum(vy)]]

# Subplot 3: Core
ax3 = axs[3] # Bottom-Left
if !isempty(core)
   vx_core = core[idx_vx, :]
   vy_core = core[idx_vy, :]
   h3 = ax3.hist2d(vx_core, vy_core, bins = 100, range = v_range,
      norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1), cmap = "viridis")
   plt.colorbar(h3[4], ax = ax3, label = "Count")
else
   ax3.text(0.5, 0.5, "No Core Particles", ha = "center", transform = ax3.transAxes)
end
ax3.set_title("Core Population")
ax3.set_xlabel("\$v_x\$")
ax3.set_ylabel("\$v_y\$")
ax3.axhline(0, color = "gray", linestyle = "--")

# Subplot 4: Suprathermal
ax4 = axs[4] # Bottom-Right
if !isempty(supra)
   vx_supra = supra[idx_vx, :]
   vy_supra = supra[idx_vy, :]
   h4 = ax4.hist2d(vx_supra, vy_supra, bins = 100, range = v_range,
      norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1), cmap = "magma")
   plt.colorbar(h4[4], ax = ax4, label = "Count")
else
   ax4.text(0.5, 0.5, "No Suprathermal Particles", ha = "center", transform = ax4.transAxes)
end
ax4.set_title("Suprathermal Population")
ax4.set_xlabel("\$v_x\$")
ax4.set_ylabel("\$v_y\$")
ax4.axhline(0, color = "gray", linestyle = "--")

# Save Unified Figure
output_filename = "phase_space_analysis.png"
savefig(output_filename)
println("Unified analysis plot saved to $output_filename")
