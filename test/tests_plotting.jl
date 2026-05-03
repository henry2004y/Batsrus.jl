using StaticArrays
@testset "Plotting" begin
    @testset "Plots" begin
        RecipesBase.is_key_supported(k::Symbol) = true
        # 1D
        file = "1d__raw_2_t25.60000_n00000258.out"
        bd = load(joinpath(datapath, file))
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bd, "Rho")
        @test getfield(rec[1], 1)[:seriestype] == :path

        file = "z=0_raw_1_t25.60000_n00000258.out"
        bd = load(joinpath(datapath, file))
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bd, "p")
        @test getfield(rec[1], 1)[:seriestype] == :contourf

        @testset "rbody and inner masking" begin
            # Mock BatsHead
            # param names: ["rbody"]
            # eqpar values: [2.5f0]
            head = Batsrus.BatsHead(
                2, "headline", 0, 0.0f0, false, 1, 1, [10, 10],
                [2.5f0], ["x", "y"], ["rho"], ["rbody"]
            )
            list = Batsrus.FileList("test", Batsrus.Real4Bat, ".", 0, 1, 0)

            # 10x10 grid from -4.5 to 4.5
            # x[i, j, 1] = -5.5 + i
            # i=5 -> -0.5, i=6 -> 0.5, i=7 -> 1.5, i=8 -> 2.5
            x = zeros(Float32, 10, 10, 2)
            for j in 1:10, i in 1:10
                x[i, j, 1] = -5.5 + i
                x[i, j, 2] = -5.5 + j
            end
            w = ones(Float32, 10, 10, 1)

            bd = Batsrus.BATS(head, list, x, w)

            @testset "Keyword overrides header" begin
                xi, yi, Wi = interp2d(bd, "rho"; innermask = true, rbody = 1.0)
                @test isnan(Wi[5, 5])
                @test isnan(Wi[6, 6])
                @test !isnan(Wi[7, 7])
            end

            @testset "Header fallback" begin
                xi, yi, Wi = interp2d(bd, "rho"; innermask = true)
                @test isnan(Wi[7, 7])
                @test !isnan(Wi[8, 8])
            end

            @testset "Default fallback" begin
                head_no_rbody = Batsrus.BatsHead(
                    2, "headline", 0, 0.0f0, false, 0, 1, [10, 10],
                    Float32[], ["x", "y"], ["rho"], String[]
                )
                bd_no_rbody = Batsrus.BATS(head_no_rbody, list, x, w)
                xi, yi, Wi = interp2d(bd_no_rbody, "rho"; innermask = true)
                @test isnan(Wi[5, 5])
                @test !isnan(Wi[7, 7])
            end
        end
    end

    if RUN_MAKIE_TESTS
        @testset "Makie" begin
            file = "1d__raw_2_t25.60000_n00000258.out"
            bd = load(joinpath(datapath, file))
            fig, ax, plt = CairoMakie.lines(bd, "Rho")
            @test plt isa Lines

            file = "z=0_raw_1_t25.60000_n00000258.out"
            bd = load(joinpath(datapath, file))
            fig, ax, plt = CairoMakie.heatmap(bd, "p")
            @test plt isa Heatmap

            @testset "Makie Animation" begin
                mktempdir() do tmpdir
                    # 2D structured binary with streamlines
                    file = "z=0_raw_1_t25.60000_n00000258.out"
                    animate(
                        [joinpath(datapath, file)],
                        outdir = tmpdir,
                        streamvars = "bx;by",
                        stream_kwargs = (; color = :red),
                        colormap = :viridis,
                        title = "Custom Title"
                    )
                    @test isfile(joinpath(tmpdir, "z=0_raw_1_t25.60000_n00000258.png"))
                end

                mktempdir() do tmpdir
                    # 2D unstructured data with streamlines
                    file_unstructured = "bx0_mhd_6_t00000100_n00000352.out"
                    animate(
                        [joinpath(datapath, file_unstructured)],
                        outdir = tmpdir,
                        var = "P", streamvars = "ux;uy",
                        stream_kwargs = (; color = :white, linewidth = 0.5)
                    )
                    @test isfile(joinpath(tmpdir, "bx0_mhd_6_t00000100_n00000352.png"))
                end
            end

            @testset "Amrex Particles" begin
                tmpdir = mktempdir()
                try
                    generate_mock_amrex_data(
                        tmpdir; num_particles = 100,
                        real_component_names = ["ux", "uy"]
                    )
                    data = AMReXParticle(tmpdir)

                    # Test plot_phase into existing axis
                    fig = CairoMakie.Figure()
                    ax = CairoMakie.Axis(fig[1, 1])
                    pl = plot_phase!(ax, data, "ux", "uy"; log_scale = false)
                    @test pl isa Makie.Plot

                    # Test plot_phase creating new figure
                    obj = plot_phase(data, "ux", "uy")
                    obj = plot_phase(data, "ux", "uy"; log_scale = false)
                    @test obj isa Makie.FigureAxisPlot
                finally
                    rm(tmpdir, recursive = true, force = true)
                end
            end
        end
    elseif RUN_PYPLOT_TESTS
        @testset "PyPlot" begin
            @test size(squeeze(zeros(2, 3, 1))) == (2, 3)
            # 1D ascii
            file = "1d__raw_2_t25.60000_n00000258.out"
            bd = load(joinpath(datapath, file), verbose = false)
            c = PyPlot.plot(bd, "p")
            @test c[1].get_xdata() ≈ bd.x
            @test c[1].get_ydata() ≈ bd.w[:, 10]

            @testset "Amrex Particles" begin
                tmpdir = mktempdir()
                try
                    generate_mock_amrex_data(
                        tmpdir; num_particles = 100,
                        real_component_names = ["ux", "uy"]
                    )
                    data = AMReXParticle(tmpdir)

                    # Test plot_phase (PyPlot creates new figure/axis by default if not provided)
                    obj = plot_phase(data, "ux", "uy")
                    # PyPlot's imshow returns a matplotlib.image.AxesImage
                    @test obj isa PyPlot.PyObject

                    # Verify we can plot into an existing axis
                    fig, ax = plt.subplots()
                    obj2 = plot_phase!(ax, data, "ux", "uy")
                    @test obj2 isa PyPlot.PyObject

                    @test ax.get_xlabel() == "ux" || ax.get_xlabel() == "\$v_x\$"
                    @test ax.get_ylabel() == "uy" || ax.get_ylabel() == "\$v_y\$"
                    plt.close(fig)

                finally
                    rm(tmpdir, recursive = true, force = true)
                end
            end

            # 2D structured binary
            file = "z=0_raw_1_t25.60000_n00000258.out"
            bd = load(joinpath(datapath, file))
            c = PyPlot.streamplot(bd, "bx;by")
            @test c.lines.get_segments()[2][3] ≈ -118.68871477694084
            c = PyPlot.contourf(bd, "p")
            @test c.get_array()[end] == 0.9750000000000002
            c = @suppress_err PyPlot.contourf(bd, "rho", innermask = true)
            @test c.get_array()[end] == 0.9750000000000002
            c = PyPlot.contour(bd, "rho")
            @test c.get_array()[end] == 1.0500000000000003
            c = PyPlot.contour(bd, "rho"; levels = [1.0])
            @test c.get_array()[end] == 1.0
            c = PyPlot.tricontourf(bd, "rho")
            @test c.get_array()[end] == 0.9750000000000002
            c = PyPlot.tricontour(bd, "rho")
            @test c isa PyPlot.PyObject
            PyPlot.tripcolor(bd, "rho")
            @test isa(gca(), PyPlot.PyObject)
            p = PyPlot.pcolormesh(bd, "p").get_array()
            @test p[end] == 0.1f0
            p = PyPlot.imshow(bd, "p").get_array()
            @test p[2, 128] == 0.51229393f0

            @testset "set_colorbar" begin
                ext = Base.get_extension(Batsrus, :BatsrusPyPlotExt)

                # twoslope
                norm_two = ext.set_colorbar(:twoslope, -1.0, 1.0, [-1.0, 1.0]; vcenter = 0.5)
                @test norm_two.vcenter == 0.5
                @test norm_two.vmin == -1.0
                @test norm_two.vmax == 1.0

                # linear
                norm_lin = ext.set_colorbar(:linear, 0.0, 1.0, [0.0, 1.0])
                @test norm_lin.vmin == 0.0
                @test norm_lin.vmax == 1.0

                # symlog
                norm_sym = ext.set_colorbar(:symlog, -10.0, 10.0, [-10.0, 10.0])
                @test norm_sym.vmin == -10.0
                @test norm_sym.vmax == 10.0

                # log
                norm_log = ext.set_colorbar(:log, 0.1, 10.0, [0.1, 10.0])
                @test norm_log.vmin == 0.1
                @test norm_log.vmax == 10.0
            end

            @testset "Streamplot Uniformity" begin
                let bd = begin
                        # Mock BatsHead
                        head = Batsrus.BatsHead(
                            2, "PLANETARY", 216, 3.0020797f0, false, 0, 2, [100, 100],
                            [], ["x", "y"], ["uxs1", "uzs1"], []
                        )

                        # Mock FileList
                        list = Batsrus.FileList(
                            "mock.out", Batsrus.Real4Bat, ".",
                            1000, 1, 0
                        )

                        # Large Float32 range that caused issues
                        x = zeros(Float32, 100, 100, 2)
                        w = rand(Float32, 100, 100, 2)

                        x[1, 1, 1] = 6.0f6
                        x[100, 1, 1] = 6.0f6 + 1.0f0
                        x[1, 1, 2] = 0.0f0
                        x[1, 100, 2] = 1.0f0

                        Batsrus.BATS(head, list, x, w)
                    end
                    # This should not throw ValueError
                    plt.figure()
                    c = PyPlot.streamplot(bd, "uxs1;uzs1")
                    @test c isa PyPlot.PyObject
                    plt.close()
                end
            end

            plt.close()
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
            c = plot_surface(bd, "rho")
            @test isa(gca(), PyPlot.PyObject)
            plt.close()

            # plot_trisurf with explicit ax
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
            c = plot_trisurf(bd, "rho", ax; cmap = "viridis")
            @test c isa PyPlot.PyObject
            plt.close()

            mktempdir() do tmpdir
                # 2D structured binary with streamlines and derived variable
                animate(
                    [joinpath(datapath, file)], outdir = tmpdir,
                    var = :b,
                    streamvars = "bx;by",
                    stream_kwargs = (; color = "red"),
                    plot_kwargs = (; cmap = "viridis"),
                    fig_kwargs = (; figsize = (4, 3)),
                    title = "Custom Title"
                )
                @test isfile(joinpath(tmpdir, "z=0_raw_1_t25.60000_n00000258.png"))
            end

            mktempdir() do tmpdir
                # 2D unstructured data with streamlines
                file_unstructured = "bx0_mhd_6_t00000100_n00000352.out"
                animate(
                    [joinpath(datapath, file_unstructured)], outdir = tmpdir,
                    var = "P", streamvars = "ux;uy",
                    stream_kwargs = (; color = "white", density = 0.5)
                )
                @test isfile(joinpath(tmpdir, "bx0_mhd_6_t00000100_n00000352.png"))
            end

            # 2D AMR Cartesian
            file = "bx0_mhd_6_t00000100_n00000352.out"
            bd = load(joinpath(datapath, file))
            pcolormesh(bd, "P")
            @test isa(gca(), PyPlot.PyObject)

            @testset "Grid Plotting" begin
                # 1. Structured BatsrusIDL
                @testset "Structured BatsrusIDL" begin
                    head = Batsrus.BatsHead(
                        2, "headline", 0, 0.0f0, false, 0, 1, [10, 10],
                        Float32[], ["x", "y"], ["rho"], String[]
                    )
                    list = Batsrus.FileList("test", Batsrus.Real4Bat, ".", 0, 1, 0)
                    x = zeros(Float32, 10, 10, 2)
                    for j in 1:10, i in 1:10
                        x[i, j, 1] = i
                        x[i, j, 2] = j
                    end
                    w = ones(Float32, 10, 10, 1)
                    bd = Batsrus.BATS(head, list, x, w)

                    plt.figure()
                    c = plotgrid(bd)
                    @test c isa PyPlot.PyObject
                    plt.close()
                end

                # 2. Batl (AMR)
                @testset "Batl (AMR)" begin
                    # Mock Head
                    head = Batsrus.Head(
                        Int32(8), Int32(8), Int32(1), # nI, nJ, nK
                        Int32(2), # nG
                        Int32(2), Int32(2), Int32(1), # iRatio, jRatio, kRatio
                        Int32(2), # nDimAmr
                        Int32(4), # nChild
                        MVector{3, Int32}(1, 1, 1), # nRoot_D
                        MVector{3, Float64}(-1.0, -1.0, -1.0), # CoordMin_D
                        MVector{3, Float64}(1.0, 1.0, 1.0), # CoordMax_D
                        Int32[1, 2], # iDimAmr_D
                        MVector{3, Bool}(false, false, false), # isPeriodic_D
                        MVector{3, Float64}(0.0, 0.0, 0.0) # dxPlot_D
                    )

                    # Mock iTree_IA
                    # 1 node (root), status=used, level=0, coord1=1, coord2=1
                    iTree_IA = zeros(Int32, 18, 1)
                    iTree_IA[Batsrus.status_, 1] = Batsrus.used_
                    iTree_IA[Batsrus.level_, 1] = 0
                    iTree_IA[Batsrus.coord1_, 1] = 1
                    iTree_IA[Batsrus.coord2_, 1] = 1

                    batl = Batl(head, iTree_IA, MVector{3, Int32}(2, 2, 1), Int8(2))

                    plt.figure()
                    plotgrid(batl)
                    @test true # Just check if it runs without error
                    plt.close()
                end

                @testset "Batl 3D (AMR)" begin
                    # Mock Head 3D
                    head = Batsrus.Head(
                        Int32(8), Int32(8), Int32(8), # nI, nJ, nK
                        Int32(2), # nG
                        Int32(2), Int32(2), Int32(2), # iRatio, jRatio, kRatio
                        Int32(3), # nDimAmr
                        Int32(8), # nChild
                        MVector{3, Int32}(1, 1, 1), # nRoot_D
                        MVector{3, Float64}(-1.0, -1.0, -1.0), # CoordMin_D
                        MVector{3, Float64}(1.0, 1.0, 1.0), # CoordMax_D
                        Int32[1, 2, 3], # iDimAmr_D
                        MVector{3, Bool}(false, false, false), # isPeriodic_D
                        MVector{3, Float64}(0.0, 0.0, 0.0) # dxPlot_D
                    )

                    # Mock iTree_IA
                    iTree_IA = zeros(Int32, 18, 1)
                    iTree_IA[Batsrus.status_, 1] = Batsrus.used_
                    iTree_IA[Batsrus.level_, 1] = 0
                    iTree_IA[Batsrus.coord1_, 1] = 1
                    iTree_IA[Batsrus.coord2_, 1] = 1
                    iTree_IA[Batsrus.coord3_, 1] = 1

                    batl = Batl(head, iTree_IA, MVector{3, Int32}(2, 2, 2), Int8(3))

                    # This will create a 3D plot
                    c = plotgrid(batl)
                    @test (c isa PyPlot.PyObject) || (c isa Vector && c[1] isa PyPlot.PyObject)
                    plt.close()

                    # 3D slice
                    plt.figure()
                    c = plotgrid(batl; dir = "x", at = 0.0)
                    @test c isa PyPlot.PyObject
                    plt.close()
                end

                # 3. Unstructured Tecplot
                @testset "Unstructured Tecplot" begin
                    head = (
                        nDim = 2,
                        nNode = 4,
                        nCell = 1,
                        variable = ["x", "y", "rho"],
                        ET = "QUADRILATERAL",
                    )
                    data = Float32[
                        0.0 1.0 1.0 0.0;
                        0.0 0.0 1.0 1.0;
                        1.0 1.0 1.0 1.0
                    ]
                    connectivity = Int32[1, 2, 3, 4][:, :]

                    plt.figure()
                    plotgrid(head, data, connectivity)
                    @test true # Just check if it runs without error
                    plt.close()
                end
            end
        end
    end
end
