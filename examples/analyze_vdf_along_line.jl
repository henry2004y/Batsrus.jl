# Script to analyze VDFs along a spatial line
# Usage: julia examples/analyze_vdf_along_line.jl

using Batsrus
using PyPlot
using Printf
using Statistics
using LinearAlgebra
using StaticArrays
using FHist

# --- Configuration ---
# Note: Update this path to your specific environment
dir_prefix = "data_mock_line_demo" # Use a separate mock data dir
data_path = dir_prefix

# Analysis Line Geometry
start_point = [17.0, -5.0, 0.0]
end_point = [19.0, -5.0, 0.0]
num_samples = 4
box_size = 0.5 # Total width of the sampling box
B_field = [1.0, -1.0, 0.0]
v_para_range = (-10.0, 10.0)
v_perp_range = (0.0, 10.0)

# ---------------------
# Section 1: Data Preparation
# ---------------------
println("--- Section 1: Data Preparation ---")

if !isdir(data_path)
   println("Data not found at: $data_path")
   println("Generating mock AMReX data for demonstration...")
   mkpath(data_path)

   # B_hat for generating data
   B_hat = normalize(B_field)

   Batsrus.generate_mock_amrex_data(
      data_path,
      num_particles = 5000000,
      real_component_names = ["ux", "uy", "uz", "weight"],
      domain_min = [-20.0, -20.0, -20.0],
      domain_max = [40.0, 40.0, 40.0],
      particle_gen = (i, n_reals) -> begin
         # Distribute particles roughly in the region of interest
         # Cylinder along X [10, 30]
         x = 15.0 + rand() * 10.0
         y = -10.0 + rand() * 10.0
         z = -5.0 + rand() * 10.0

         # Spatially varying distribution
         # Near start (17.0): Core dominated
         # Near end (19.0): Beam dominated
         eff_pos = (x - 17.0) / 2.0 # 0 to 1 roughly
         beam_prob = clamp(eff_pos, 0.1, 0.9) # simple linear mix

         is_beam = rand() < beam_prob

         if is_beam
            # Fast Beam aligned with B
            v_para = 5.0 + randn() * 1.0
            v_perp1 = randn() * 2.0
            v_perp2 = randn() * 2.0
         else
            # Static Core
            v_para = randn() * 1.5
            v_perp1 = randn() * 1.5
            v_perp2 = randn() * 1.5
         end

         # Convert "local aligned" v to global Cartesian
         # Construct basis
         b = B_hat
         # Random perp for basis
         rand_vec = [1.0, 0.0, 0.0]
         if abs(dot(rand_vec, b)) > 0.9
            rand_vec = [0.0, 1.0, 0.0]
         end
         perp1 = normalize(cross(b, rand_vec))
         perp2 = cross(b, perp1)

         v_vec = (v_para .* b) .+ (v_perp1 .* perp1) .+ (v_perp2 .* perp2)

         return (
            x, y, z,
            v_vec[1], v_vec[2], v_vec[3],
            1.0
         )
      end
   )
   println("Mock data generated at: $data_path")
end

# Load data
data = AMReXParticle(data_path)
println("Data loaded. Dimensions: $(data.dim)D")

# ---------------------
# Section 2: Analysis Loop
# ---------------------
println("\n--- Section 2: Extracting VDFs along the line ---")

# Calculate sampling points
deltas = (end_point .- start_point) ./ (num_samples - 1)
sample_points = [start_point .+ (i - 1) .* deltas for i in 1:num_samples]

# Prepare plot
fig, axs = plt.subplots(
   1, num_samples, figsize = (4 * num_samples, 4), constrained_layout = true)

# Transformation function
transform_func = get_particle_field_aligned_transform(B_field)

plot_data = []
global_max_density = 1.0

println("Pass 1: Identifying global density range...")
for (i, pt) in enumerate(sample_points)
   println("  Sample $i: $pt")

   # Define box around point
   hw = box_size / 2
   x_range = (pt[1] - hw, pt[1] + hw)
   y_range = (pt[2] - hw, pt[2] + hw)
   z_range = (pt[3] - hw, pt[3] + hw)

   # 2D/3D selection support
   local particles
   if data.dim == 2
      particles = select_particles_in_region(data; x_range, y_range)
   else
      particles = select_particles_in_region(data; x_range, y_range, z_range)
   end

   n_part = size(particles, 2)

   if n_part < 10
      push!(plot_data, nothing)
      continue
   end

   t_data, t_names = transform_func(particles, data.header.real_component_names)
   v_para = t_data[1, :]
   v_perp = t_data[2, :]

   # Pre-calculate histogram to find max density
   h = Hist2D((v_para, v_perp),
      binedges = (range(v_para_range..., length = 51), range(v_perp_range..., length = 51)))
   counts = bincounts(h)
   local_max = maximum(counts)

   global global_max_density = max(global_max_density, local_max)

   push!(plot_data, counts)
end

println("Global Max Density: $global_max_density")

println("Pass 2: Plotting...")
for (i, pt) in enumerate(sample_points)
   ax = axs[i]
   data_pair = plot_data[i]

   if isnothing(data_pair)
      ax.text(0.5, 0.5, "Not enough particles", ha = "center", va = "center")
      ax.set_title("Point $i")
      continue
   end

   counts = data_pair

   # Plot with fixed global normalization
   # Transpose counts because imshow expects (row, col) -> (y, x)
   # but bincounts returns (x, y)
   im = ax.imshow(counts',
      extent = [v_para_range[1], v_para_range[2], v_perp_range[1], v_perp_range[2]],
      origin = "lower",
      norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1, vmax = global_max_density),
      cmap = "turbo",
      aspect = "auto",
      interpolation = "nearest")

   ax.set_title(@sprintf("Pt %d: x=%.1f", i, pt[1]))
   ax.set_xlabel(L"v_{\parallel}")
   if i == 1
      ax.set_ylabel(L"v_{\perp}")
   end

   # Save last handle for colorbar
   global h_plot = im
end

if @isdefined(h_plot)
   # Create a dummy ScalarMappable with the global norm for the colorbar
   norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1, vmax = global_max_density)
   sm = PyPlot.matplotlib.cm.ScalarMappable(norm = norm, cmap = "turbo")
   sm.set_array([])
   fig.colorbar(sm, ax = axs, label = "Counts", pad = 0.02)
end

# Add a global title
fig.suptitle("VDF Evolution Along Line: $(start_point) -> $(end_point)")

outfile = "vdf_line_analysis.png"
savefig(outfile)
println("\nAnalysis complete. Plot saved to $outfile")
