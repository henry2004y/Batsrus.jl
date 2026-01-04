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
box_size = 1.0 # Total width of the sampling box
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
         beam_prob = clamp(eff_pos * 0.5, 0.0, 0.4) # limit beam to 40% to keep core density dominant

         is_beam = rand() < beam_prob

         if is_beam
            # Fast Beam aligned with B
            v_para = 5.0 + randn() * 1.0
            v_perp1 = randn() * 2.0
            v_perp2 = randn() * 2.0
         else
            # Static Core
            # Make core slightly colder (sigma 1.2) to ensure high phase space density
            v_para = randn() * 1.2
            v_perp1 = randn() * 1.2
            v_perp2 = randn() * 1.2
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

# Transformation function
base_transform = get_particle_field_aligned_transform(B_field)

transform_func = (data, names) -> begin
   t_data, t_names = base_transform(data, names)

   # Extract parallel velocity and weights
   # v_para is row 1, v_perp is row 2, weight is row 3
   v_para = @view t_data[1, :]
   weights = @view t_data[3, :]

   # Find bulk parallel velocity (peak location)
   if !isempty(v_para)
      vmin, vmax = extrema(v_para)
      if vmin < vmax
         # 100 bins for estimation is sufficient
         h1 = Hist1D(v_para; binedges = range(vmin, vmax, length = 101), weights = weights)
         _, max_idx = findmax(h1.bincounts)
         edges = h1.binedges isa Tuple ? h1.binedges[1] : h1.binedges
         bulk_v_para = (edges[max_idx] + edges[max_idx + 1]) / 2

         # Shift parallel velocity to rest frame
         v_para .-= bulk_v_para
      end
   end

   return t_data, t_names
end

plot_data = []
global_max_density = 1e-25 # Small initial value for PSD

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

   # Calculate histogram using get_phase_space_density
   # density calculation is now default
   edges = (range(v_para_range..., length = 51), range(v_perp_range..., length = 51))
   h = get_phase_space_density(
      data, "v_parallel", "v_perp";
      x_range, y_range, z_range,
      edges,
      transform = transform_func
   )
   psd_code = Float64.(h.bincounts)

   # Convert from code units to s^3 km^-6
   # Code density is N / (L_code^3 * V_code^3)  (assuming values are km/s already)
   # We need N / (L_phys^3 * V_phys^3)
   # L_phys = L_code * RE_km
   # V_phys = V_code
   # So psd_phys = psd_code / RE_km^3

   RE_km = 6378.0 # km

   psd = psd_code ./ (RE_km^3)

   local_max = maximum(psd)

   global global_max_density = max(global_max_density, local_max)

   push!(plot_data, psd)
end

println("Global Max Density: $global_max_density")

println("Pass 2: Plotting...")

fig, axs = plt.subplots(
   1, num_samples, figsize = (4 * num_samples, 4), constrained_layout = true)

for (i, pt) in enumerate(sample_points)
   ax = axs[i]
   psd = plot_data[i]

   if isnothing(psd)
      ax.text(0.5, 0.5, "Not enough particles", ha = "center", va = "center")
      ax.set_title("Point $i")
      continue
   end

   # Plot PSD
   # Transpose psd because imshow expects (row, col) -> (y, x)
   # Determine vmin dynamically to handle potentially small SI values
   vmin_val = max(1e-40, global_max_density * 1e-6)

   im = ax.imshow(psd',
      extent = [v_para_range[1], v_para_range[2], v_perp_range[1], v_perp_range[2]],
      origin = "lower",
      norm = PyPlot.matplotlib.colors.LogNorm(vmin = vmin_val, vmax = global_max_density),
      cmap = "turbo",
      aspect = "auto",
      interpolation = "nearest")

   ax.set_title(@sprintf("Pt %d: x=%.1f", i, pt[1]))
   ax.set_xlabel(L"v_{\parallel} \ [km/s]")
   if i == 1
      ax.set_ylabel(L"v_{\perp} \ [km/s]")
   else
      ax.set_yticklabels([])
   end

   # Save last handle for colorbar
   global h_plot = im
end

if @isdefined(h_plot)
   # Create a dummy ScalarMappable with the global norm for the colorbar
   vmin_val = max(1e-40, global_max_density * 1e-6)
   norm = PyPlot.matplotlib.colors.LogNorm(vmin = vmin_val, vmax = global_max_density)
   sm = PyPlot.matplotlib.cm.ScalarMappable(norm = norm, cmap = "turbo")
   sm.set_array([])
   fig.colorbar(sm, ax = axs, label = L"Phase Space Density $[s^3 km^{-6}]$", pad = 0.02)
end

# Add a global title
fig.suptitle("VDF Evolution Along Line: $(start_point) -> $(end_point)")

outfile = "vdf_line_analysis.png"
savefig(outfile)
println("\nAnalysis complete. Plot saved to $outfile")
