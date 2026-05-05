module test_show

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
        end
    end
end

end # module test_show
