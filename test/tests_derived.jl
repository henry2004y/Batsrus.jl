@testset "Derived Variables and Analysis" begin
    # 2D structured binary
    file_struct = "z=0_raw_1_t25.60000_n00000258.out"
    bd_struct = load(joinpath(datapath, file_struct))

    @testset "Field Magnitudes and Vectors" begin
        # B field
        @test Batsrus.fill_vector_from_scalars(bd_struct, :B)[:, end, end] ==
            Float32[1.118034, -0.559017, 0.0]

        b_mag = get_magnitude(bd_struct, :B)
        @test b_mag[128, 2] == 0.9223745f0
        @test bd_struct[:b] == b_mag

        b_mag2 = get_magnitude2(bd_struct, :B)
        @test b_mag2[128, 2] == 0.8507747f0
        @test bd_struct[:b2] == b_mag2
    end

    # 2D region/unstructured binary
    file_aniso = "z=0_fluid_region0_0_t00001640_n00010142.out"
    bd_aniso = load(joinpath(datapath, file_aniso))

    @testset "Electric and Velocity Fields" begin
        # E field
        e_mag = get_magnitude(bd_aniso, :E)
        @test e_mag[2, 1] == 2655.4805f0
        @test bd_aniso[:e] == e_mag

        @test get_magnitude2(bd_aniso, :E)[2, 1] == 7.051577f6
        @test Batsrus.fill_vector_from_scalars(bd_aniso, :E)[:, 2, 1] ==
            Float32[-241.05942, -2644.2058, -40.53219]

        # U field - this file uses species suffixes (uxs0, etc.)
        @test get_magnitude2(bd_aniso, :U0)[2, 1] == 33784.973f0
        # For this file, :u (which calls :U) will fail because 'ux' doesn't exist.
        # We test get_magnitude with :U0 instead.
        u0_mag = get_magnitude(bd_aniso, :U0)
        @test u0_mag[2, 1] == √33784.973f0
    end

    @testset "Anisotropy" begin
        a0 = get_anisotropy(bd_aniso, 0)
        @test a0[1:2, 1] ≈ Float32[1.2630985, 2.4700143]
        @test bd_aniso[:anisotropy0] == a0

        # Test rotation method
        @test get_anisotropy(bd_aniso, 0, method = :rotation)[1:2, 1] ≈ a0[1:2, 1]

        a1 = get_anisotropy(bd_aniso, 1)
        @test a1[1:2, 1] ≈ Float32[1.2906302, 2.6070855]
        @test bd_aniso[:anisotropy1] == a1
    end

    @testset "Extended Electrodynamics" begin
        w_conv = get_convection_E(bd_aniso)
        @test w_conv[2][2, 1] ≈ -2454.3933f0

        w_hall = get_hall_E(bd_aniso)
        @test w_hall[2][2, 1] ≈ -782.2945f0
    end

    @testset "Time Series Analysis" begin
        w = get_timeseries([joinpath(datapath, file_aniso)], [0.0, 0.0])
        @test w[2][end] == 17.973747f0
    end

    @testset "Unstructured ASCII Derived" begin
        file_ascii = "bx0_mhd_6_t00000100_n00000352.out"
        bd_ascii = load(joinpath(datapath, file_ascii))
        @test get_magnitude(bd_ascii, :U)[2] == 71.85452748407637
        @test get_magnitude2(bd_ascii, :U)[2] == 5163.073119959886
        @test bd_ascii[:u] == get_magnitude(bd_ascii, :U)
    end

    @testset "Fallback behavior" begin
        @test bd_struct[:rho] == bd_struct["rho"]
    end
end
