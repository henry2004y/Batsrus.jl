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
        @testset "Structured 2D Mock" begin
            nx, ny = 2, 2
            names = [
                SubString("bx"), SubString("by"), SubString("bz"),
                SubString("pxxs0"), SubString("pyys0"), SubString("pzzs0"),
                SubString("pxys0"), SubString("pxzs0"), SubString("pyzs0"),
                SubString("pxxs1"), SubString("pyys1"), SubString("pzzs1"),
                SubString("pxys1"), SubString("pxzs1"), SubString("pyzs1"),
            ]
            w_data = zeros(Float32, nx, ny, length(names))

            # Magnetic field along X: B = [1, 0, 0]
            w_data[:, :, 1] .= 1.0f0

            # Species 0: P = [2 0 0; 0 1 0; 0 0 1] -> p_parallel = 2, p_perp = 1 -> A = 0.5
            w_data[:, :, 4] .= 2.0f0 # pxx
            w_data[:, :, 5] .= 1.0f0 # pyy
            w_data[:, :, 6] .= 1.0f0 # pzz

            # Species 1: P = [1 0 0; 0 2 0; 0 0 2] -> p_parallel = 1, p_perp = 2 -> A = 2.0
            w_data[:, :, 10] .= 1.0f0 # pxx
            w_data[:, :, 11] .= 2.0f0 # pyy
            w_data[:, :, 12] .= 2.0f0 # pzz

            head = Batsrus.BatsHead(
                2, SubString("Mock"), 0, 0.0f0, false, 0, length(names), [nx, ny], Float32[],
                [SubString("x"), SubString("y")], names, String[]
            )
            x_data = zeros(Float32, nx, ny, 2)
            bd = BATS(
                head,
                Batsrus.FileList("mock_aniso.out", Batsrus.Real4Bat, ".", 0, 1, 0),
                x_data, w_data
            )

            # Test Projection Method
            a0 = get_anisotropy(bd, 0, method = :projection)
            a1 = get_anisotropy(bd, 1, method = :projection)
            @test all(a0 .== 0.5f0)
            @test all(a1 .== 2.0f0)

            # Test Rotation Method
            ar0 = get_anisotropy(bd, 0, method = :rotation)
            ar1 = get_anisotropy(bd, 1, method = :rotation)
            @test all(isapprox.(ar0, 0.5f0, atol = 1.0e-6))
            @test all(isapprox.(ar1, 2.0f0, atol = 1.0e-6))

            # Test Symbol Access
            @test all(isapprox.(bd[:anisotropy0], 0.5f0, atol = 1.0e-6))
            @test all(isapprox.(bd[:anisotropy1], 2.0f0, atol = 1.0e-6))
        end

        @testset "Structured 3D Rotation Mock" begin
            # Test rotation with non-trivial B field
            nx, ny, nz = 2, 2, 2
            names = [
                SubString("bx"), SubString("by"), SubString("bz"),
                SubString("pxxs0"), SubString("pyys0"), SubString("pzzs0"),
                SubString("pxys0"), SubString("pxzs0"), SubString("pyzs0"),
            ]
            w_data = zeros(Float32, nx, ny, nz, length(names))

            # B = [1, 1, 1]
            w_data[:, :, :, 1] .= 1.0f0
            w_data[:, :, :, 2] .= 1.0f0
            w_data[:, :, :, 3] .= 1.0f0

            # P = I (identity) -> A = 1.0
            w_data[:, :, :, 4] .= 1.0f0 # pxx
            w_data[:, :, :, 5] .= 1.0f0 # pyy
            w_data[:, :, :, 6] .= 1.0f0 # pzz

            head = Batsrus.BatsHead(
                3, SubString("Mock"), 0, 0.0f0, false, 0, length(names), [nx, ny, nz], Float32[],
                [SubString("x"), SubString("y"), SubString("z")], names, String[]
            )
            x_data = zeros(Float32, nx, ny, nz, 3)
            bd = BATS(
                head,
                Batsrus.FileList("mock_aniso3d.out", Batsrus.Real4Bat, ".", 0, 1, 0),
                x_data, w_data
            )

            ar = get_anisotropy(bd, 0, method = :rotation)
            @test all(isapprox.(ar, 1.0f0, atol = 1.0e-6))
        end

        @testset "Curvilinear Grid Mock" begin
            # For ndim=1, BATS expects (nx, 1) for x_data
            nx = 4
            names = [
                SubString("bx"), SubString("by"), SubString("bz"),
                SubString("pxxs0"), SubString("pyys0"), SubString("pzzs0"),
                SubString("pxys0"), SubString("pxzs0"), SubString("pyzs0"),
            ]
            w_data = zeros(Float32, nx, length(names))
            w_data[:, 1] .= 1.0f0 # bx = 1
            w_data[:, 4] .= 2.0f0 # pxx = 2
            w_data[:, 5] .= 1.0f0 # pyy = 1
            w_data[:, 6] .= 1.0f0 # pzz = 1

            head = Batsrus.BatsHead(
                1, SubString("Mock"), 0, 0.0f0, true, 0, length(names), [nx], Float32[],
                [SubString("x")], names, String[]
            )
            x_data = zeros(Float32, nx, 2)
            bd = BATS(
                head,
                Batsrus.FileList("mock_aniso_u.out", Batsrus.Real4Bat, ".", 0, 1, 0),
                x_data, w_data
            )

            a0 = get_anisotropy(bd, 0, method = :projection)
            @test all(a0 .== 0.5f0)

            ar0 = get_anisotropy(bd, 0, method = :rotation)
            @test all(isapprox.(ar0, 0.5f0, atol = 1.0e-6))
        end

        @testset "Real Data" begin
            file_aniso = "z=0_fluid_region0_0_t00001640_n00010142.out"
            bd_aniso = load(joinpath(datapath, file_aniso))

            a0 = get_anisotropy(bd_aniso, 0)
            @test a0[1:2, 1] ≈ Float32[1.2630985, 2.4700143]
            @test bd_aniso[:anisotropy0] == a0

            # Test rotation method
            @test get_anisotropy(bd_aniso, 0, method = :rotation)[1:2, 1] ≈ a0[1:2, 1]

            a1 = get_anisotropy(bd_aniso, 1)
            @test a1[1:2, 1] ≈ Float32[1.2906302, 2.6070855]
            @test bd_aniso[:anisotropy1] == a1
        end
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

    @testset "Current Density" begin
        # 2D structured verification
        nx, ny = 5, 5
        names = ["bx", "by", "bz"]
        nw = length(names)
        w_data = zeros(Float32, nx, ny, nw)
        # B_y = x, so J_z = dBy/dx = 1
        x_range = range(0.0f0, 4.0f0, length = nx) # dx = 1.0
        y_range = range(0.0f0, 4.0f0, length = ny) # dy = 1.0
        for iy in 1:ny, ix in 1:nx
            w_data[ix, iy, 2] = x_range[ix]
        end

        head = Batsrus.BatsHead(
            2, "Mock", 0, 0.0, false, 0, nw, [nx, ny], Float32[],
            ["x", "y"], names, String[]
        )
        x_data = zeros(Float32, nx, ny, 2)
        for iy in 1:ny, ix in 1:nx
            x_data[ix, iy, 1] = x_range[ix]
            x_data[ix, iy, 2] = y_range[iy]
        end
        list = Batsrus.FileList("mock_j.out", Batsrus.Real4Bat, ".", 0, 1, 0)
        bd_mock = BATS(head, list, x_data, w_data)

        jx, jy, jz = get_current_density(bd_mock)
        @test all(jx .== 0.0f0)
        @test all(jy .== 0.0f0)
        @test all(jz .≈ 1.0f0)

        @test bd_mock[:jz] == jz
        @test bd_mock[:j] ≈ sqrt.(jx .^ 2 .+ jy .^ 2 .+ jz .^ 2)
    end
end
