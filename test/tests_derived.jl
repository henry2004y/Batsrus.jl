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

        # 3D Structured Mock Magnitude
        nx, ny, nz = 2, 2, 2
        w_data = zeros(Float32, nx, ny, nz, 3)
        w_data[:, :, :, 1] .= 1.0f0 # bx
        w_data[:, :, :, 2] .= 1.0f0 # by
        w_data[:, :, :, 3] .= 1.0f0 # bz
        head = Batsrus.BatsHead(
            3, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx, ny, nz], Float32[],
            [SubString("x"), SubString("y"), SubString("z")],
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        bd3d = BATS(head, Batsrus.FileList("mock_3d.out", Batsrus.Real4Bat, ".", 0, 1, 0), zeros(Float32, nx, ny, nz, 3), w_data)
        @test all(get_magnitude(bd3d, :B) .≈ √3.0f0)

        # 4D Generic Mock Magnitude (Fallback case)
        sz = [2, 2, 2, 2]
        w_data_nd = zeros(Float32, sz..., 3)
        w_data_nd[:, :, :, :, 1] .= 1.0f0
        head_nd = Batsrus.BatsHead(
            4, SubString("Mock"), 0, 0.0f0, false, 0, 3, sz, Float32[],
            [SubString("x1"), SubString("x2"), SubString("x3"), SubString("x4")],
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        bd_nd = BATS(head_nd, Batsrus.FileList("mock_4d.out", Batsrus.Real4Bat, ".", 0, 1, 0), zeros(Float32, sz..., 4), w_data_nd)
        @test all(get_magnitude(bd_nd, :B) .== 1.0f0)
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

            # Test species 2 (tests _get_pressure_tensor_indices species 2 branch)
            # We mock a file with species 0, 1, 2
            names_s2 = [
                "bx", "by", "bz",
                "pxxs0", "pyys0", "pzzs0", "pxys0", "pxzs0", "pyzs0",
                "pxxs1", "pyys1", "pzzs1", "pxys1", "pxzs1", "pyzs1",
                "pxxs2", "pyys2", "pzzs2", "pxys2", "pxzs2", "pyzs2",
            ]
            w_data_s2 = zeros(Float32, 2, length(names_s2))
            w_data_s2[:, 1] .= 1.0f0 # bx=1
            w_data_s2[:, 16] .= 3.0f0 # pxxs2=3
            w_data_s2[:, 17] .= 1.0f0 # pyys2=1
            w_data_s2[:, 18] .= 1.0f0 # pzzs2=1

            head_s2 = Batsrus.BatsHead(
                1, "Mock", 0, 0.0, true, 0, length(names_s2), [2], Float32[],
                ["x"], names_s2, String[]
            )
            bd_s2 = BATS(head_s2, Batsrus.FileList("mock_s2.out", Batsrus.Real4Bat, ".", 0, 1, 0), zeros(Float32, 2, 2), w_data_s2)
            @test all(get_anisotropy(bd_s2, 2) .≈ 1 / 3.0f0)

            # Test species 3 (tests generic species string branch)
            names_s3 = [
                "bx", "by", "bz",
                "pxxs0", "pyys0", "pzzs0", "pxys0", "pxzs0", "pyzs0",
                "pxxs1", "pyys1", "pzzs1", "pxys1", "pxzs1", "pyzs1",
                "pxxs2", "pyys2", "pzzs2", "pxys2", "pxzs2", "pyzs2",
                "pxxs3", "pyys3", "pzzs3", "pxys3", "pxzs3", "pyzs3",
            ]
            w_data_s3 = zeros(Float32, 2, length(names_s3))
            w_data_s3[:, 1] .= 1.0f0
            w_data_s3[:, 22] .= 4.0f0 # pxxs3=4
            w_data_s3[:, 23] .= 1.0f0 # pyys3=1
            w_data_s3[:, 24] .= 1.0f0 # pzzs3=1
            head_s3 = Batsrus.BatsHead(
                1, "Mock", 0, 0.0, true, 0, length(names_s3), [2], Float32[],
                ["x"], names_s3, String[]
            )
            bd_s3 = BATS(head_s3, Batsrus.FileList("mock_s3.out", Batsrus.Real4Bat, ".", 0, 1, 0), zeros(Float32, 2, 2), w_data_s3)
            @test all(get_anisotropy(bd_s3, 3) .≈ 0.25f0)
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
        @testset "1D Structured" begin
            # B_y = x, so J_z = dBy/dx = 1
            nx = 5
            x_range = range(0.0f0, 4.0f0, length = nx) # dx = 1.0
            w_data = zeros(Float32, nx, 3)
            w_data[:, 2] = x_range # by

            head = Batsrus.BatsHead(
                1, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx], Float32[],
                [SubString("x")], [SubString("bx"), SubString("by"), SubString("bz")], String[]
            )
            x_data = collect(reshape(x_range, nx, 1))
            list = Batsrus.FileList("mock_1d.out", Batsrus.Real4Bat, ".", 0, 1, 0)
            bd = BATS(head, list, x_data, w_data)

            jx, jy, jz = get_current_density(bd)
            @test all(jx .== 0.0f0)
            @test all(jy .== -0.0f0)
            @test jz[2:4] ≈ [1.0f0, 1.0f0, 1.0f0]

            # Test individual component access (triggers _compute_j*)
            @test all(bd[:jx] .== 0.0f0)
            @test all(bd[:jy] .== 0.0f0)
            @test bd[:jz][2:4] ≈ [1.0f0, 1.0f0, 1.0f0]
        end

        @testset "2D Structured" begin
            # B_z = x, so J_y = -dBz/dx = -1
            # B_z = y, so J_x = dBz/dy = 1
            nx, ny = 5, 5
            x_range = range(0.0f0, 4.0f0, length = nx)
            y_range = range(0.0f0, 4.0f0, length = ny)

            w_data = zeros(Float32, nx, ny, 3)
            for iy in 1:ny, ix in 1:nx
                w_data[ix, iy, 3] = x_range[ix] + y_range[iy] # bz = x + y
            end

            head = Batsrus.BatsHead(
                2, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx, ny], Float32[],
                [SubString("x"), SubString("y")],
                [SubString("bx"), SubString("by"), SubString("bz")], String[]
            )
            # BATS expects (nx, ny, 2) for x_data in 2D
            x_data = zeros(Float32, nx, ny, 2)
            for iy in 1:ny, ix in 1:nx
                x_data[ix, iy, 1] = x_range[ix]
                x_data[ix, iy, 2] = y_range[iy]
            end
            bd = BATS(head, Batsrus.FileList("mock_2d.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

            jx, jy, jz = get_current_density(bd)
            @test jx[2:4, 2:4] ≈ ones(Float32, 3, 3)    # dBz/dy = 1
            @test jy[2:4, 2:4] ≈ -ones(Float32, 3, 3)   # -dBz/dx = -1
            @test all(jz .== 0.0f0)

            # Test individual component access
            @test bd[:jx][2:4, 2:4] ≈ ones(Float32, 3, 3)
            @test bd[:jy][2:4, 2:4] ≈ -ones(Float32, 3, 3)
            @test all(bd[:jz] .== 0.0f0)

            # Test total current density
            @test bd[:j] ≈ sqrt.(bd[:jx] .^ 2 .+ bd[:jy] .^ 2 .+ bd[:jz] .^ 2)
        end

        @testset "3D Structured" begin
            nx, ny, nz = 5, 5, 5
            x_range = range(0.0f0, 4.0f0, length = nx)
            y_range = range(0.0f0, 4.0f0, length = ny)
            z_range = range(0.0f0, 4.0f0, length = nz)

            w_data = zeros(Float32, nx, ny, nz, 3)
            for iz in 1:nz, iy in 1:ny, ix in 1:nx
                w_data[ix, iy, iz, 1] = y_range[iy] # bx = y
                w_data[ix, iy, iz, 2] = x_range[ix] # by = x
            end

            head = Batsrus.BatsHead(
                3, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx, ny, nz], Float32[],
                [SubString("x"), SubString("y"), SubString("z")],
                [SubString("bx"), SubString("by"), SubString("bz")], String[]
            )
            x_data = zeros(Float32, nx, ny, nz, 3)
            for iz in 1:nz, iy in 1:ny, ix in 1:nx
                x_data[ix, iy, iz, 1] = x_range[ix]
                x_data[ix, iy, iz, 2] = y_range[iy]
                x_data[ix, iy, iz, 3] = z_range[iz]
            end
            bd = BATS(head, Batsrus.FileList("mock_3d.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

            jx, jy, jz = get_current_density(bd)
            @test all(jx .== 0.0f0)
            @test all(jy .== 0.0f0)
            @test all(isapprox.(jz[2:4, 2:4, 2:4], 0.0f0, atol = 1.0e-6))

            # Test individual component access
            @test all(bd[:jx] .== 0.0f0)
            @test all(bd[:jy] .== 0.0f0)
            @test all(isapprox.(bd[:jz][2:4, 2:4, 2:4], 0.0f0, atol = 1.0e-6))
        end

        @testset "Unstructured/Curvilinear Grid Error" begin
            # For ndim=1, BATS expects (nx, 1) for x_data
            nx = 4
            head = Batsrus.BatsHead(
                1, SubString("Mock"), 0, 0.0f0, true, 0, 3, [nx], Float32[],
                [SubString("x")],
                [SubString("bx"), SubString("by"), SubString("bz")], String[]
            )
            x_data = zeros(Float32, nx, 2) # Dimp1 = 2 = ndim+1
            w_data = zeros(Float32, nx, 3)
            bd = BATS(head, Batsrus.FileList("mock_u.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

            @test_throws Exception get_current_density(bd)
        end

        @testset "PLANETARY Units Scaling" begin
            nx, ny = 5, 5
            x_range = range(0.0f0, 4.0f0, length = nx)
            y_range = range(0.0f0, 4.0f0, length = ny)
            w_data = zeros(Float32, nx, ny, 3)
            for iy in 1:ny, ix in 1:nx
                w_data[ix, iy, 2] = x_range[ix]
            end

            head = Batsrus.BatsHead(
                2, SubString("PLANETARY"), 0, 0.0f0, false, 0, 3, [nx, ny], Float32[],
                [SubString("x"), SubString("y")],
                [SubString("bx"), SubString("by"), SubString("bz")], String[]
            )
            x_data = zeros(Float32, nx, ny, 2)
            for iy in 1:ny, ix in 1:nx
                x_data[ix, iy, 1] = x_range[ix]
                x_data[ix, iy, 2] = y_range[iy]
            end
            bd = BATS(head, Batsrus.FileList("mock_p.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

            jx, jy, jz = get_current_density(bd)
            @test jz[2, 2] ≈ Float32(Batsrus.FAC_J_PLANETARY)
        end

        @testset "Existing J in file" begin
            # Test _getvar(:j) when jx, jy, jz already exist
            nx, ny = 2, 2
            names = ["jx", "jy", "jz"]
            w_data = zeros(Float32, nx, ny, 3)
            w_data[:, :, 1] .= 1.0f0 # jx
            w_data[:, :, 2] .= 1.0f0 # jy
            w_data[:, :, 3] .= 1.0f0 # jz
            head = Batsrus.BatsHead(
                2, "Mock", 0, 0.0, false, 0, 3, [nx, ny], Float32[],
                ["x", "y"], names, String[]
            )
            bd = BATS(head, Batsrus.FileList("mock_j_exist.out", Batsrus.Real4Bat, ".", 0, 1, 0), zeros(Float32, nx, ny, 2), w_data)
            @test all(bd[:j] .≈ √3.0f0)
        end
    end

    @testset "Electron Pressure Gradient E-field" begin
        # 2D structured mock data
        nx, ny = 5, 5
        x_range = range(0.0f0, 1.0f0, length = nx)
        y_range = range(0.0f0, 1.0f0, length = ny)
        # rhos0, pxxs0, pyys0, pzzs0, pxys0, pxzs0, pyzs0
        w_data = zeros(Float32, nx, ny, 7)
        w_data[:, :, 1] .= 1.0f0 # rhos0
        for ix in 1:nx, iy in 1:ny
            w_data[ix, iy, 2] = x_range[ix] # pxxs0 = x -> dPxx/dx = 1
        end
        head = Batsrus.BatsHead(
            2, "normalized", 0, 0.0f0, false, 0, 7, [nx, ny], Float32[],
            ["x", "y"],
            ["rhos0", "pxxs0", "pyys0", "pzzs0", "pxys0", "pxzs0", "pyzs0"], String[]
        )
        x_data = zeros(Float32, nx, ny, 2)
        for ix in 1:nx, iy in 1:ny
            x_data[ix, iy, 1] = x_range[ix]
            x_data[ix, iy, 2] = y_range[iy]
        end
        bd = BATS(head, Batsrus.FileList("mock_pe.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

        Ex, Ey, Ez = get_pe_E(bd)
        @test size(Ex) == (nx, ny)
        # C=1, ne=1, dPxx/dx=1, other divP components 0
        @test all(Ex .≈ -1.0f0)
        @test all(Ey .≈ 0.0f0)
        @test all(Ez .≈ 0.0f0)

        # Test manual mass override: ne = rho / m. If m=2, ne=0.5, Ex=-2
        Ex_m, _, _ = get_pe_E(bd, mass = 2.0)
        @test all(Ex_m .≈ -2.0f0)

        # 3D structured mock data
        nx, ny, nz = 4, 4, 4
        x_range = range(0.0f0, 1.0f0, length = nx)
        y_range = range(0.0f0, 1.0f0, length = ny)
        z_range = range(0.0f0, 1.0f0, length = nz)
        w_data3d = zeros(Float32, nx, ny, nz, 7)
        w_data3d[:, :, :, 1] .= 1.0f0 # rhos0
        for ix in 1:nx, iy in 1:ny, iz in 1:nz
            w_data3d[ix, iy, iz, 2] = x_range[ix] # pxxs0 = x
        end
        head3d = Batsrus.BatsHead(
            3, "normalized", 0, 0.0f0, false, 0, 7, [nx, ny, nz], Float32[],
            ["x", "y", "z"],
            ["rhos0", "pxxs0", "pyys0", "pzzs0", "pxys0", "pxzs0", "pyzs0"], String[]
        )
        x_data3d = zeros(Float32, nx, ny, nz, 3)
        for ix in 1:nx, iy in 1:ny, iz in 1:nz
            x_data3d[ix, iy, iz, 1] = x_range[ix]
            x_data3d[ix, iy, iz, 2] = y_range[iy]
            x_data3d[ix, iy, iz, 3] = z_range[iz]
        end
        bd3d = BATS(head3d, Batsrus.FileList("mock_pe3d.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data3d, w_data3d)
        Ex3d, Ey3d, Ez3d = get_pe_E(bd3d)
        @test size(Ex3d) == (nx, ny, nz)
        @test all(Ex3d .≈ -1.0f0)
    end
end
