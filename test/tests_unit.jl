@testset "Units" begin
    # Mockup for PLANETARY units
    # Using NamedTuple to avoid struct definition issues in test scope
    head_mock = (
        headline = "PLANETARY", wname = ["x", "rho", "ux", "bx", "p", "jx"], ndim = 3,
    )
    bd_mock = (head = head_mock,)

    @test 1.0u"Re" > 2.0u"Rg"
    @test 1.0u"RSun" > 1.0u"Re"
    @test 1.0u"Re" > 1.0u"RMercury"
    @test 1.0u"AU" > 1.0u"RSun"
    file = "y=0_var_1_t00000000_n00000000.out"
    bd = load(joinpath(datapath, file))
    varunit = getunit(bd, "Rho")
    @test varunit == u"amucc"
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

    # Mockup for PLANETARY units

    @test getunit(bd_mock, "x") == u"Re"
    @test getunit(bd_mock, "rho") == u"amucc"
    @test getunit(bd_mock, "ux") == u"km/s"
    @test getunit(bd_mock, "bx") == u"nT"
    @test getunit(bd_mock, "p") == u"nPa"
    @test getunit(bd_mock, "jx") == u"μA/m^2"

    # Test getunits
    units = getunits(bd_mock)
    @test length(units) == 6
    @test units[1] == u"Re"
    @test units[2] == u"amucc"
    @test units[6] == u"μA/m^2"
end
