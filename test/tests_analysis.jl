using Random

@testset "Particle Classification" begin
    mktempdir() do tmpdir
        Random.seed!(1234)
        n_core = 1000
        n_halo = 100
        n_total = n_core + n_halo

        Batsrus.generate_mock_amrex_data(
            tmpdir;
            num_particles = n_total,
            real_component_names = ["vx", "vy", "vz"],
            particle_gen = (i, n_reals) -> begin
                # Positions (random in box 0-2)
                pos = rand(3) .* 2.0 .- 0.5

                # Velocities
                if i <= n_core
                    # Core velocities (Gaussian)
                    vel = randn(3)
                else
                    # Halo velocities (Broad Gaussian)
                    vel = randn(3) .* 5.0
                end

                (pos..., vel...)
            end
        )

        data = AMReXParticle(tmpdir)
        # Ensure data is loaded
        Batsrus.load_data!(data)

        # Test 1: Classification with known params
        # vth=1.0, u=0
        # nsigma=3.0

        core, halo = classify_particles(
            data; vdim = 3, bulk_vel = [0.0, 0.0, 0.0], vth = 1.0, nsigma = 3.0
        )

        # Just check counts roughly
        @test size(core, 2) > n_core * 0.9
        @test size(halo, 2) < n_halo + n_core * 0.1
        @test size(halo, 2) > n_halo * 0.5 # Halo shouldn't be empty

        # Test 2: Auto validation of bulk velocity
        core_auto, halo_auto = classify_particles(
            data; vdim = 3, bulk_vel = nothing, vth = 1.0, nsigma = 3.0
        )

        @test size(core_auto, 2) ≈ size(core, 2) atol = 50

        # Test 3: 1D
        core_1d, halo_1d = classify_particles(
            data; vdim = 1, bulk_vel = [0.0], vth = 1.0, nsigma = 3.0
        )

        # In 1D, we verify against x-velocity only.
        # Similar statistics hold.
        @test size(core_1d, 2) > n_core * 0.9

        # Test 4: Direct test of get_core_population_mask
        velocities = [
            0.0 1.0 2.0 10.0;
            0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0
        ]
        bulk_vel = [0.0, 0.0, 0.0]
        vth = 1.0
        nsigma = 3.0
        mask = Batsrus.get_core_population_mask(velocities, bulk_vel, vth, nsigma)

        @test mask == [true, true, true, false]
    end
end

@testset "Field Aligned Transform" begin
    # Test with simple magnetic field along X
    # v_para should be vx, v_perp should be sqrt(vy^2 + vz^2)
    b_field = [1.0, 0.0, 0.0]
    transform = Batsrus.get_particle_field_aligned_transform(b_field)

    # Mock data: 2 particles, 6 components (x, y, z, vx, vy, vz)
    # Particle 1: vx=1, vy=0, vz=0 => v_para=1, v_perp=0
    # Particle 2: vx=0, vy=3, vz=4 => v_para=0, v_perp=5
    names = ["x", "y", "z", "vx", "vy", "vz"]
    data = [
        0.0 0.0 0.0 1.0 0.0 0.0;
        0.0 0.0 0.0 0.0 3.0 4.0
    ]'

    new_data, new_names = transform(data, names)

    @test new_names == ["v_parallel", "v_perp", "weight"]
    @test new_data[1, 1] ≈ 1.0
    @test new_data[1, 2] ≈ 0.0
    @test new_data[2, 1] ≈ 0.0
    @test new_data[2, 2] ≈ 5.0

    # Test with E and B
    # B = x, E = y
    # ExB = z
    # b_hat = (1,0,0)
    # d_hat = (0,0,1)
    # e_hat = d x b = (0,0,1) x (1,0,0) = (0,1,0) (which is E direction)
    # So basis is (x, y, z) corresponding to (v_B, v_E, v_BxE)

    b_field = [10.0, 0.0, 0.0]
    e_field = [0.0, 5.0, 0.0]
    transform_eb = Batsrus.get_particle_field_aligned_transform(b_field, e_field)

    # Particle: vx=1, vy=2, vz=3
    # v_B = 1
    # v_E = 2
    # v_BxE = 3
    data_eb = [0.0 0.0 0.0 1.0 2.0 3.0]'

    new_data_eb, new_names_eb = transform_eb(data_eb, names)

    @test new_names_eb == ["v_B", "v_E", "v_BxE", "weight"]
    @test new_data_eb[1, 1] ≈ 1.0
    @test new_data_eb[2, 1] ≈ 2.0
    @test new_data_eb[3, 1] ≈ 3.0
end

@testset "GMM Fit" begin
    mktempdir() do tmpdir
        Random.seed!(5678)
        # Generate mock data: 2 populations
        # Population 1 (Core): Weight 0.7, Mean 0, vth 1.0
        # Population 2 (Beam): Weight 0.3, Mean 5.0, vth 0.5

        n1 = 7000
        n2 = 3000
        n_total = n1 + n2

        # Use Float32 to test generic support
        T = Float32

        s1 = T(1.0 / sqrt(2))
        s2 = T(0.5 / sqrt(2))

        Batsrus.generate_mock_amrex_data(
            tmpdir;
            num_particles = n_total,
            real_component_names = ["vx", "vy", "vz"],
            particle_gen = (i, n_reals) -> begin
                pos = rand(3)

                if i <= n1
                    vel = randn(T, 3) .* s1
                else
                    vel = randn(T, 3) .* s2
                    vel[1] += 5.0f0 # Shift vx
                end

                (pos..., vel...)
            end
        )

        data = AMReXParticle(tmpdir)

        # Fit 2 clusters
        results = @suppress fit_particle_velocity_gmm(data, 2)

        @test length(results) == 2

        # Check parameters (Results sorted by weight)
        c1 = results[1] # Expected Core
        @test isapprox(c1.weight, 0.7f0, atol = 0.05)
        @test isapprox(c1.mean[1], 0.0f0, atol = 0.2)
        @test isapprox(c1.vth[1], 1.0f0, atol = 0.2)

        c2 = results[2] # Expected Beam
        @test isapprox(c2.weight, 0.3f0, atol = 0.05)
        @test isapprox(c2.mean[1], 5.0f0, atol = 0.2)
        @test isapprox(c2.vth[1], 0.5f0, atol = 0.2)

        # Test Thermal Velocity Extraction
        b_field = [1.0f0, 0.0f0, 0.0f0]
        # c1 (core) is isotropic with vth=1.0. Parallel should be ~1.0, perp ~1.0.
        vpara, vperp = Batsrus.get_gmm_thermal_velocity(c1, b_field)
        @test isapprox(vpara, 1.0f0, atol = 0.2)
        @test isapprox(vperp, 1.0f0, atol = 0.2)
    end
end
