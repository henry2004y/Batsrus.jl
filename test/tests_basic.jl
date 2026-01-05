
@testset "Units" begin
   @test 1.0bu"R" > 2.0bu"Rg"
   file = "y=0_var_1_t00000000_n00000000.out"
   bd = load(joinpath(datapath, file))
   varunit = getunit(bd, "Rho")
   @test varunit == bu"amucc"
   varunit = getunit(bd, "Ux")
   @test varunit == u"km/s"
   varunit = getunit(bd, "Bx")
   @test varunit == u"nT"
   varunit = getunit(bd, "P")
   @test varunit == u"nPa"
   varunit = getunit(bd, "jx")
   @test dimension(varunit) == dimension(Unitful.A / Unitful.m^2)
   varunit = getunit(bd, "ex")
   @test varunit == u"mV/m"
end
