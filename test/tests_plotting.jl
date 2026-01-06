
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
   end

   @testset "Makie" begin
      file = "1d__raw_2_t25.60000_n00000258.out"
      bd = load(joinpath(datapath, file))
      fig, ax, plt = CairoMakie.lines(bd, "Rho")
      @test plt isa Lines

      file = "z=0_raw_1_t25.60000_n00000258.out"
      bd = load(joinpath(datapath, file))
      fig, ax, plt = CairoMakie.heatmap(bd, "p")
      @test plt isa Heatmap

      @testset "Amrex Particles" begin
         tmpdir = mktempdir()
         try
            generate_mock_amrex_data(tmpdir; num_particles = 100,
               real_component_names = ["ux", "uy"])
            data = AMReXParticle(tmpdir)

            # Test plot_phase into existing axis
            fig = CairoMakie.Figure()
            ax = CairoMakie.Axis(fig[1, 1])
            pl = plot_phase(data, "ux", "uy"; ax = ax)
            @test pl isa Makie.Plot

            # Test plot_phase creating new figure
            obj = plot_phase(data, "ux", "uy")
            @test obj isa Makie.FigureAxisPlot
         finally
            rm(tmpdir, recursive = true, force = true)
         end
      end
   end

   if RUN_PYPLOT_TESTS
      @testset "PyPlot" begin
         @test size(squeeze(zeros(2, 3, 1))) == (2, 3)
         # 1D ascii
         file = "1d__raw_2_t25.60000_n00000258.out"
         bd = load(joinpath(datapath, file), verbose = false)
         c = PyPlot.plot(bd, "p")
         @test c[1].get_xdata() ≈ bd.x
         @test c[1].get_ydata() ≈ bd.w[:, 10]

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
         PyPlot.tripcolor(bd, "rho")
         @test isa(gca(), PyPlot.PyObject)
         p = PyPlot.pcolormesh(bd, "p").get_array()
         @test p[end] == 0.1f0
         p = PyPlot.imshow(bd, "p").get_array()
         @test p[2, 128] == 0.51229393f0
         plt.close()
         fig = plt.figure()
         ax = fig.add_subplot(111, projection = "3d")
         plot_surface(bd, "rho")
         @test isa(gca(), PyPlot.PyObject)
         plt.close()

         # 2D AMR Cartesian
         file = "bx0_mhd_6_t00000100_n00000352.out"
         bd = load(joinpath(datapath, file))
         pcolormesh(bd, "P")
         @test isa(gca(), PyPlot.PyObject)
      end
   end
end
