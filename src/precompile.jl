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
         Batsrus.generate_mock_amrex_data(tmpdir)
         data = AMReXParticle(tmpdir)
         get_phase_space_density(data, "x", "u"; bins = 2)
      end
   end
end
