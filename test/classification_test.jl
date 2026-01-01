
using Batsrus, Test
using FHist

@testset "Particle Classification" begin
   mktempdir() do tmpdir
      Batsrus.generate_mock_amrex_data(tmpdir)
      data = AMReXParticle(tmpdir)

      # Ensure data is loaded to populate _rdata
      Batsrus.load_data!(data)

      # Modify header to include velocity aliases if not present
      # Mock data has minimal columns. We will OVERWRITE _rdata with 6 columns.
      # And we must trick the system to think these columns exist.
      # data.header.real_component_names is a Vector{String}.

      # Check what's in there currently
      # println("Current components: ", data.header.real_component_names)

      # Resize component names to 6 if needed
      # We just replace it.
      empty!(data.header.real_component_names)
      append!(data.header.real_component_names, ["x", "y", "z", "vx", "vy", "vz"])

      n_core = 1000
      n_halo = 100
      n_total = n_core + n_halo

      # Create synthetic data
      rdata = zeros(Float64, n_total, 6)

      # Positions (random in box 0-2 (since mock data is small?))
      # Mock data domain dimensions: [2, 2, 2].
      # Left edge: -0.5, Right edge: 1.5. Length = 2.0.
      rdata[:, 1:3] .= rand(n_total, 3) .* 2.0 .- 0.5

      # Core velocities (Gaussian)
      rdata[1:n_core, 4:6] .= randn(n_core, 3)

      # Halo velocities (Broad Gaussian)
      rdata[(n_core + 1):end, 4:6] .= randn(n_halo, 3) .* 5.0

      data._rdata = rdata

      # Test 1: Classification with known params
      # vth=1.0, u=0
      # nsigma=3.0

      core, halo = classify_particles(
         data; vdim = 3, bulk_vel = [0.0, 0.0, 0.0], vth = 1.0, nsigma = 3.0)

      # Core efficiency: erte of 3 sigma is 0.997.
      # So expected core count ~ n_core * 0.997 + halo_leak
      # Halo leak: halo sigma=5. threshold=3. fraction within 3/5=0.6 sigma of halo distribution?
      # erf(0.6/sqrt(2)) is small.

      # Just check counts roughly
      @test size(core, 1) > n_core * 0.9
      @test size(halo, 1) < n_halo + n_core * 0.1
      @test size(halo, 1) > n_halo * 0.5 # Halo shouldn't be empty

      # Test 2: Auto validation of bulk velocity
      core_auto, halo_auto = classify_particles(
         data; vdim = 3, bulk_vel = nothing, vth = 1.0, nsigma = 3.0)

      @test size(core_auto, 1)≈size(core, 1) atol=50

      # Test 3: 1D
      core_1d, halo_1d = classify_particles(
         data; vdim = 1, bulk_vel = [0.0], vth = 1.0, nsigma = 3.0)

      # In 1D, we verify against x-velocity only.
      # Similar statistics hold.
      @test size(core_1d, 1) > n_core * 0.9
   end
end
