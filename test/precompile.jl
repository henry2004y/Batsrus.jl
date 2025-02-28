# Precompiling workloads

@setup_workload begin
   @compile_workload begin
      # Get a minimal file for precompilation!
      file = "z=0_raw_1_t25.60000_n00000258.out"
      bd = load(joinpath(datapath, file))
      plotrange = [-10.0, 10.0, -Inf, Inf]
      x, y, w = interp2d(bd, "rho", plotrange)
      Batsrus.fill_vector_from_scalars(bd, :B)
      get_magnitude(bd, :B)
      get_magnitude2(bd, :B)
      # Linear interpolation at a given point
      d = interp1d(bd, "rho", Float32[0.0, 0.0])
      # Linear interpolation along a line
      point1 = Float32[-10.0, -1.0]
      point2 = Float32[10.0, 1.0]
      w = interp1d(bd, "rho", point1, point2)
      w = slice1d(bd, "rho", 1, 1)
      get_var_range(bd, "rho")
   end
end