using FortranFiles

@testset "Reading 1D ascii" begin
    file = "1d__raw_2_t25.60000_n00000258.out"
    bd = @suppress_err begin
        load(joinpath(datapath, file), verbose = true)
    end
    @test startswith(repr(bd), "filename : 1d")
    @test extrema(bd.x) == (-127.5, 127.5)
    @test extrema(bd.w) == (-0.79960780498, 1.9394335293)
end

@testset "Reading 2D structured binary" begin
    file = "z=0_raw_1_t25.60000_n00000258.out"
    @test_throws ArgumentError load(joinpath(datapath, file), npict = 2)
    bd = load(joinpath(datapath, file))
    @test bd.head.time == 25.6f0
    @test extrema(bd.x) == (-127.5f0, 127.5f0)
    @test extrema(bd.w) == (-0.79985905f0, 1.9399388f0)
    plotrange = [-10.0, 10.0, -Inf, Inf]
    x, y, w = interp2d(bd, "rho", plotrange)
    @test w[1, end] ≈ 0.6848635077476501 atol = 1.0e-4
    # Linear interpolation at a given point
    d = interp1d(bd, "rho", Float32[0.0, 0.0])
    @test d == 0.6936918f0
    # Linear interpolation along a line
    point1 = Float32[-10.0, -1.0]
    point2 = Float32[10.0, 1.0]
    w = interp1d(bd, "rho", point1, point2)
    @test sum(w) == 14.942542f0
    w = slice1d(bd, "rho", 1, 1)
    @test sum(w) == 4.0f0
    @test get_var_range(bd, "rho") == (0.11626893f0, 1.0f0)
    @test length(bd["rho"][X = -10 .. 10, Y = Near(0.9)]) == 20
    @test bd["rho"][X = At(0.5), Y = At(0.5)] == 0.6940014f0

    file = "z=0_fluid_region0_0_t00001640_n00010142.out"
    bd = load(joinpath(datapath, file))
    x, y = Batsrus.meshgrid(bd)
    @test length(x) == 601 && y[2] == 0.0f0
    x, y = Batsrus.meshgrid(bd, Float32[-100, 100, -Inf, Inf])
    @test length(x) == 4
end

@testset "Reading 2D unstructured ascii" begin
    file = "bx0_mhd_6_t00000100_n00000352.out"
    bd = load(joinpath(datapath, file))
    plotrange = [-Inf, Inf, -Inf, Inf]
    x, y = Batsrus.meshgrid(bd)
    @test length(x) == 117 && length(y) == 246
    x, y, w = interp2d(bd, "rho", plotrange, useMatplotlib = false)
    @test w[1, 2] == 5.000018304080387
end

@testset "Reading 3D structured binary" begin
    file = "3d_raw.out"
    bd = load(joinpath(datapath, file))
    plotrange = [-50.0, 50.0, -0.5, 0.5]
    X, Z, p = cutdata(bd, "p"; dir = "y", sequence = 1, plotrange)
    @test p[1] ≈ 0.560976f0 && p[2] ≈ 0.53704995f0
    @test size(bd["p"]) == (8, 8, 8)
end

@testset "Log" begin
    logfilename = joinpath(datapath, "log_n000001.log")
    head, bd = readlogdata(logfilename)
    @test isa(head, NamedTuple)
    @test extrema(bd) == (-0.105, 258.0)
end

@testset "VTK" begin
    using WriteVTK, Glob, LightXML
    file = joinpath(datapath, "z=0_fluid_region0_0_t00001640_n00010142.out")
    convertIDLtoVTK(file)
    sha_str = bytes2hex(open(sha1, "out.vti"))
    @test sha_str == "7c7540e684aaf1823b7e2792f555e89ad98bcd31"

    file = joinpath(datapath, "3d_raw.out")
    convertIDLtoVTK(file)
    sha_str = bytes2hex(open(sha1, "out.vti"))
    @test sha_str == "8d535c0656daa9224aab048d036da0cc2f667bdc"
    rm("out.vti")

    file = joinpath(datapath, "3d_bin.dat")
    convertTECtoVTU(file)
    sha_str = bytes2hex(open(sha1, "out.vtu"))
    @test sha_str == "5b04747666542d802357dec183177f757754a254"
    rm("out.vtu")

    filetag = joinpath(datapath, "3d_mhd_amr/3d__mhd_1_t00000000_n00000000")
    batl = Batl(readhead(filetag * ".info"), readtree(filetag)...)
    # local block index check
    @test Batsrus.find_grid_block(batl, [1.0, 0.0, 0.0]) == 4 &&
        Batsrus.find_grid_block(batl, [100.0, 0.0, 0.0]) == -100
    connectivity = getConnectivity(batl)
    sha_str = bytes2hex(sha256(string(connectivity)))
    @test sha_str == "c6c5a65a46d86a9ba4096228c1516f89275e45e295cd305eb70c281a770ede74"
end

@testset "HDF5" begin
    using HDF5
    file = "3d__var_3_n00000000_single.batl"
    bd = BatsrusHDF5Uniform(joinpath(datapath, file))
    # non-lazy
    var, _, _ = extract_var(bd, "bx")
    @test size(var) == (8, 8, 4) && eltype(var) == Float32
    # lazy
    var_lazy, _, _ = extract_var(bd, "bx"; lazy = true)
    @test var == var_lazy

    file = "3d__var_3_n00000000_double.batl"
    bd = BatsrusHDF5Uniform(joinpath(datapath, file))
    # non-lazy
    var, _, _ = extract_var(bd, "bx")
    @test size(var) == (8, 8, 4) && eltype(var) == Float64

    # Test close
    @test isopen(bd.common.fid)
    close(bd)
    @test !isopen(bd.common.fid)
end

@testset "AMReX Loader" begin
    mktempdir() do tmpdir
        Batsrus.generate_mock_amrex_data(tmpdir)

        # Test Loading
        data = AMReXParticle(tmpdir)
        @test data.dim == 3
        @test data.header.num_levels == 1
        @test data.header.num_particles == 10

        # Test Data Access
        rdata = data.rdata
        @test size(rdata) == (5, 10)
        @test rdata[1, 1] == 1.0 # x of first particle
        @test rdata[4, 1] == 10.0 # u of first particle
        @test rdata[5, 10] == 1000.0 # v of last particle

        # Test Region Selection
        selected = select_particles_in_region(data, x_range = (2.5, 4.5))
        @test size(selected, 2) == 2
        @test selected[1, 1] == 3.0
        @test selected[1, 2] == 4.0

        # Test Plotting Helper
        h = get_phase_space_density(data, "x", "u", bins = 2)
        @test size(h.bincounts) == (2, 2)

        # Test Weighted Plotting
        Batsrus.generate_mock_amrex_data(
            tmpdir;
            num_particles = 10,
            real_component_names = ["u", "v", "weight"],
            particle_gen = (i, n_reals) -> (
                Float64(i), 10.0, 10.0, # x, y, z
                Float64(i), 0.0, # u, v
                Float64(2.0), # weight
            )
        )
        data_w = AMReXParticle(tmpdir)

        # Helper to calculate bin volume
        get_bin_vol(h) = prod(e -> e[2] - e[1], h.binedges)

        h_w = get_phase_space_density(
            data_w, "x", "u", bins = 1,
            x_range = (0.0, 11.0), y_range = (0.0, 11.0)
        )

        # 10 particles, each weight 2.0. Total sum should be 20.0.
        # Density check: sum(density) * bin_vol * spatial_vol ≈ Total Weight
        # x_range=11, y_range=11, z_range default=10 (domain size)
        vol_spatial = 11.0 * 11.0 * 10.0
        bin_vol = get_bin_vol(h_w)

        @test sum(h_w.bincounts) * bin_vol * vol_spatial ≈ 20.0

        # Test 1D Plotting
        h1 = get_phase_space_density(data_w, "x"; bins = 1, x_range = (0.0, 11.0))
        vol_spatial_1d = 11.0 * 10.0 * 10.0 # x=11, y=10, z=10 implicitly
        bin_vol_1d = get_bin_vol(h1)

        @test sum(h1.bincounts) * bin_vol_1d * vol_spatial_1d ≈ 20.0
        @test ndims(h1.bincounts) == 1

        # Test 1D density
        hist1d = get_phase_space_density(data, "x")
        @test ndims(hist1d.bincounts) == 1
        # Mock data has 10 particles, weight 1 (default). Total 10.
        # Spatial vol = 10*10*10 = 1000.
        # bin_vol...

        # Test 1D density with z_range filter (Regression test)
        # Filter z to a range that includes particles (z=2.0 and z=3.0)
        # Particles are at z = i for i=1..10.
        # z_range (1.5, 3.5) captures i=2 and i=3. Total 2 particles, weight 1.
        # x_range (1.5, 3.5) captures i=2 and i=3.
        # Wait, particle i has x=i, y=i, z=i.
        # So i=2: x=2, y=2, z=2.
        # i=3: x=3, y=3, z=3.
        # Both satisfy ranges.
        # Total count = 2.

        hist1d_z = get_phase_space_density(
            data, "x";
            z_range = (1.5, 3.5), x_range = (1.5, 3.5)
        )

        # Spatial vol = (3.5-1.5) * 10 * (3.5-1.5) = 2 * 10 * 2 = 40.0
        # (y range is default 10)
        vol_spatial_z = 2.0 * 10.0 * 2.0
        bin_vol_z = get_bin_vol(hist1d_z)

        @test sum(hist1d_z.bincounts) * bin_vol_z * vol_spatial_z ≈ 2.0

        # Filter z to a range that EXCLUDES particles
        # This might error if empty, check behavior.
        # Current implementation errors on empty: "No particles found..."
        @test_throws ErrorException get_phase_space_density(
            data, "x"; z_range = (11.0, 20.0)
        )

        # Test 3D Plotting
        h3 = get_phase_space_density(
            data_w, "x", "y", "u"; bins = 1,
            x_range = (0.0, 11.0), y_range = (0.0, 20.0), z_range = (0.0, 11.0)
        )

        # Spatial vol = 11 * 20 * 11
        vol_spatial_3d = 11.0 * 20.0 * 11.0
        bin_vol_3 = get_bin_vol(h3)

        @test sum(h3.bincounts) * bin_vol_3 * vol_spatial_3d ≈ 20.0
        @test ndims(h3.bincounts) == 3
    end
end

@testset "Issue 103: readtree with header/ASCII" begin
    mktempdir() do tmpdir
        filetag = joinpath(tmpdir, "test")

        @testset "ASCII with header" begin
            open(filetag * ".tree", "w") do f
                println(f, "BATL tree information")
                println(f, "#START")
                println(f, "          3          18           2")
                println(f, "          2           2           2")
                println(f, "          1           1           1")
                for i in 1:36
                    print(f, i, " ")
                    if i % 18 == 0
                        println(f)
                    end
                end
            end
            iTree_IA, iRatio_D, nDim = readtree(filetag)
            @test nDim == 3
            @test size(iTree_IA) == (18, 2)
            @test iTree_IA[1, 1] == 1
            @test iTree_IA[18, 2] == 36
        end

        @testset "Regular binary" begin
            f = FortranFile(filetag * ".tree", "w")
            write(f, Int32(3), Int32(18), Int32(1))
            write(f, Int32[2, 2, 2])
            write(f, Int32[1, 1, 1])
            write(f, fill(Int32(7), (18, 1)))
            close(f)
            iTree_IA, iRatio_D, nDim = readtree(filetag)
            @test nDim == 3
            @test size(iTree_IA) == (18, 1)
            @test iTree_IA[1, 1] == 7
        end

        @testset "Binary with header" begin
            open(filetag * ".tree", "w") do io
                println(io, "Custom Header line 1")
                println(io, "Custom Header line 2")
                println(io, "#START")
                # Manually write Fortran records
                # Record 1: nDim, nInfo, nNode (3 * 4 = 12 bytes)
                write(io, Int32(12))
                write(io, Int32(3), Int32(18), Int32(1))
                write(io, Int32(12))
                # Record 2: iRatio_D (3 * 4 = 12 bytes)
                write(io, Int32(12))
                write(io, Int32[2, 2, 2])
                write(io, Int32(12))
                # Record 3: nRoot_D (3 * 4 = 12 bytes)
                write(io, Int32(12))
                write(io, Int32[1, 1, 1])
                write(io, Int32(12))
                # Record 4: iTree_IA (18 * 1 * 4 = 72 bytes)
                write(io, Int32(72))
                write(io, fill(Int32(9), (18, 1)))
                write(io, Int32(72))
            end
            iTree_IA, iRatio_D, nDim = readtree(filetag)
            @test nDim == 3
            @test iTree_IA[1, 1] == 9
        end

        @testset "ASCII without header" begin
            open(filetag * ".tree", "w") do f
                println(f, "          3          18           1")
                println(f, "          2           2           2")
                println(f, "          1           1           1")
                for i in 1:18
                    print(f, i, " ")
                end
                println(f)
            end
            iTree_IA, iRatio_D, nDim = readtree(filetag)
            @test nDim == 3
            @test iTree_IA[1, 1] == 1
            @test iTree_IA[18, 1] == 18
        end
    end
end
