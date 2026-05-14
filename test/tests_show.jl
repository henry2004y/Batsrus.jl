using Test
using Batsrus
using Suppressor

@testset "Show Methods" begin
    if @isdefined(datapath)
        file = joinpath(datapath, "1d__raw_2_t25.60000_n00000258.out")
        if isfile(file)
            bd = load(file)

            @testset "Compact show" begin
                s = repr(bd)
                @test startswith(s, "BatsrusIDL(1D, Float64")
            end

            @testset "Detailed show" begin
                s = repr("text/plain", bd)
                @test contains(s, "Batsrus Data Header")
                @test contains(s, "filename")
                @test contains(s, "filetype")
                @test contains(s, "iteration")
                @test contains(s, "filesize")
                @test contains(s, "snapshots")
            end

            @testset "showhead" begin
                s = @capture_out showhead(bd)
                @test contains(s, "coordinates")
                @test contains(s, "variables")
            end

            @testset "Batl show" begin
                filetag = joinpath(datapath, "3d_mhd_amr/3d__mhd_1_t00000000_n00000000")
                if isfile(filetag * ".tree")
                    batl = Batl(readhead(filetag * ".info"), readtree(filetag)...)
                    s = repr(batl)
                    @test contains(s, "Batl(3D,")
                    s_plain = repr("text/plain", batl)
                    @test contains(s_plain, "BATL 3D AMR grid")
                    @test contains(s_plain, "Nodes:")
                end
            end
        end
    end
end
