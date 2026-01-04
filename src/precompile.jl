# Precompiling workloads
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
   @compile_workload begin
      # Get a minimal file for precompilation!
      file = joinpath(@__DIR__, "../test/precompile.out")
      bd = file |> load
      x, y, w = interp2d(bd, "rho")
      Batsrus.fill_vector_from_scalars(bd, :B)
      get_magnitude(bd, :B)
      get_magnitude2(bd, :B)
      # Linear interpolation at a given point
      d = interp1d(bd, "rho", Float32[0.0, 0.0])
      # Linear interpolation along a line
      point1 = Float32[-1.0, -1.0]
      point2 = Float32[1.0, 1.0]
      w = interp1d(bd, "rho", point1, point2)
      w = slice1d(bd, "rho", 1, 1)
      get_var_range(bd, "rho")

      mktempdir() do tmpdir
         Batsrus.generate_mock_amrex_data(tmpdir;
            real_component_names = ["ux", "uy", "uz"],
            particle_gen = (i, n_reals) -> (
               rand(), rand(), rand(), randn(), randn(), randn())
         )
         data = AMReXParticle(tmpdir)
         get_phase_space_density(data, "x", "ux"; bins = 2)

         particles = select_particles_in_region(data)
         transform_func = Batsrus.get_particle_field_aligned_transform([1.0, 0.0, 0.0])
         transformed_data, _ = transform_func(
            particles, data.header.real_component_names)

         vels = transformed_data[1:3, :]
         Batsrus.get_core_population_mask(vels, [0.0, 0.0, 0.0], 1.0)

         Batsrus.fit_particle_velocity_gmm(vels, 1)
      end
   end
end
