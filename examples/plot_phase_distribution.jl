using Batsrus
using PyPlot
using Printf

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

# Load the data
data = AMReXParticle(data_path)

# Define ranges
half_width = width / 2
x_range = (x_center - half_width, x_center + half_width)
z_range = (z_center - half_width, z_center + half_width)

println("Selecting particles in region:")
println("  x: $x_range")
println("  z: $z_range")

# Plot vx - vy phase space
# Note: 'vx' and 'vy' are standard aliases for 'velocity_x' and 'velocity_y'

# Handle dimension mapping for 2D X-Z simulations
# If the data is 2D, Batsrus.jl reads it as (x, y).
# For an X-Z simulation, 'y' in the data structure usually corresponds to 'z' in physical space.

# 1. Determine ranges based on dimensions
if data.dim == 2
   println("Detected 2D data. Asuming X-Z simulation, mapping z_range to y_range.")
   # Map z constraint to 2nd dimension for selection
   selection_y_range = z_range
else
   selection_y_range = nothing # Or actual y_range if we had one
end

# 2. Select particles
# Note: For 3D, we would pass z_range. For 2D X-Z, we pass the mapped y_range.
# We are strictly following the user request to call `select_particles_in_region`
if data.dim == 2
   particles = select_particles_in_region(
      data; x_range = x_range, y_range = selection_y_range)
else
   particles = select_particles_in_region(data; x_range = x_range, z_range = z_range)
end

println("Number of particles selected: ", size(particles, 2))

# 3. Identify velocity components
names = data.header.real_component_names
# Helper to find index, handling potential aliases if needed, though here we look for standard ones
function find_component_index(names, target_aliases)
   for (i, name) in enumerate(names)
      if name in target_aliases
         return i
      end
   end
   error("Component $(target_aliases) not found in data.")
end

idx_vx = find_component_index(names, ["vx", "velocity_x", "ux"])
idx_vy = find_component_index(names, ["vy", "velocity_y", "uy"])

# Extract data
vx = particles[idx_vx, :]
vy = particles[idx_vy, :]

fig = figure(figsize = (10, 8))
# 4. Plot
h = plt.hist2d(
   vx, vy, bins = 100, norm = PyPlot.matplotlib.colors.LogNorm(), cmap = "viridis")
plt.colorbar(label = "Count")

# Add zero lines
plt.axhline(0, color = "gray", linestyle = "--", linewidth = 1, alpha = 0.5)
plt.axvline(0, color = "gray", linestyle = "--", linewidth = 1, alpha = 0.5)

xlabel("v_x")
ylabel("v_y")

title("Region: x=$(x_range), z=$(z_range)")

# Save the figure
output_filename = "phase_space_vx_vy.png"
savefig(output_filename)
println("Plot saved to $output_filename")

# --- Classification Demo ---

using Statistics

println("\n--- Classification Demo ---")

# Estimate thermal velocity from the data (using standard deviation of vx, vy, vz)
# We need to find vz index as well for 3D estimation, or just use what we have.
idx_vz = find_component_index(names, ["vz", "velocity_z", "uz"])
vz = particles[idx_vz, :]

# Calculate thermal velocity (scalar estimate)
vth_x = std(vx)
vth_y = std(vy)
vth_est = sqrt((vth_x^2 + vth_y^2) / 2)

println("Estimated thermal velocity (vth): ", vth_est)

# Classify particles
# We use the same ranges (already filtered, but the function takes ranges to filter again if needed.
# Since we already have `data` loaded, we can just pass the ranges again.
# Note: classify_particles re-selects from `data` based on ranges.
println("Classifying particles...")

core, supra = classify_particles(
   data;
   x_range = x_range,
   y_range = selection_y_range, # Handle 2D mapping
   vdim = 3,
   vth = vth_est,
   nsigma = 2.0 # Use 2 sigma for demonstration to ensure we get some suprathermals
)

println("Core particles: ", size(core, 2))
println("Suprathermal particles: ", size(supra, 2))

# Plot Classification Results
# Use subplots directly
fig, axs = plt.subplots(
   1, 2; figsize = (12, 5), constrained_layout = true, sharex = true, sharey = true)

# Calculate common range for consistent binning and plotting
v_range = [[minimum(vx), maximum(vx)], [minimum(vy), maximum(vy)]]

# Subplot 1: Core
ax1 = axs[1]
if !isempty(core)
   vx_core = core[idx_vx, :]
   vy_core = core[idx_vy, :]
   # Note: ax.hist2d returns (counts, xedges, yedges, image)
   h1 = ax1.hist2d(vx_core, vy_core, bins = 100, range = v_range,
      norm = PyPlot.matplotlib.colors.LogNorm(), cmap = "viridis")
   # h1[4] is the image object (QuadMesh), needed for colorbar
   plt.colorbar(h1[4], ax = ax1, label = "Count")
else
   ax1.text(0.5, 0.5, "No Core Particles", ha = "center", transform = ax1.transAxes)
end
ax1.set_title("Core Population")
ax1.set_xlabel("\$v_x\$")
ax1.set_ylabel("\$v_y\$")
ax1.axhline(0, color = "gray", linestyle = "--")
ax1.axvline(0, color = "gray", linestyle = "--")

# Subplot 2: Suprathermal
ax2 = axs[2]
if !isempty(supra)
   vx_supra = supra[idx_vx, :]
   vy_supra = supra[idx_vy, :]
   h2 = ax2.hist2d(vx_supra, vy_supra, bins = 100, range = v_range,
      norm = PyPlot.matplotlib.colors.LogNorm(), cmap = "magma")
   plt.colorbar(h2[4], ax = ax2, label = "Count")
else
   ax2.text(0.5, 0.5, "No Suprathermal Particles", ha = "center", transform = ax2.transAxes)
end
ax2.set_title("Suprathermal Population")
ax2.set_xlabel("\$v_x\$")
ax2.set_ylabel("\$v_y\$")
ax2.axhline(0, color = "gray", linestyle = "--")
ax2.axvline(0, color = "gray", linestyle = "--")

output_filename_class = "phase_space_classified.png"
savefig(output_filename_class)
println("Classification plot saved to $output_filename_class")
