module test_current_density

using Test
using Batsrus
using DimensionalData

@testset "Current Density" begin
    @testset "1D Structured" begin
        # B_y = x, so J_z = dBy/dx = 1
        nx = 5
        x_range = range(0.0f0, 4.0f0, length = nx) # dx = 1.0
        w_data = zeros(Float32, nx, 3)
        w_data[:, 2] = x_range # by

        # Use SubString to match BatsHead definition exactly
        head = Batsrus.BatsHead(
            1, SubString("Mock"), 0, 0.0, false, 0, 3, [nx], Float32[],
            [SubString("x")], [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        x_data = collect(reshape(x_range, nx, 1))
        list = Batsrus.FileList("mock_1d.out", Batsrus.Real4Bat, ".", 0, 1, 0)
        bd = BATS(head, list, x_data, w_data)

        jx, jy, jz = get_current_density(bd)
        @test all(jx .== 0.0f0)
        @test all(jy .== -0.0f0) # ∂Bz/∂x = 0
        @test jz[2:4] ≈ [1.0f0, 1.0f0, 1.0f0] # central diff (2.0 - 0.0)/2.0 = 1.0
    end

    @testset "3D Structured" begin
        # B_y = x, B_x = y, B_z = 0
        # J_z = dBy/dx - dBx/dy = 1 - 1 = 0
        nx, ny, nz = 5, 5, 5
        x_range = range(0.0f0, 4.0f0, length = nx) # dx = 1.0
        y_range = range(0.0f0, 4.0f0, length = ny) # dy = 1.0
        z_range = range(0.0f0, 4.0f0, length = nz) # dz = 1.0

        w_data = zeros(Float32, nx, ny, nz, 3)
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            w_data[ix, iy, iz, 1] = y_range[iy] # bx
            w_data[ix, iy, iz, 2] = x_range[ix] # by
        end

        head = Batsrus.BatsHead(
            3, SubString("Mock"), 0, 0.0, false, 0, 3, [nx, ny, nz], Float32[],
            [SubString("x"), SubString("y"), SubString("z")],
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        x_data = zeros(Float32, nx, ny, nz, 3)
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            x_data[ix, iy, iz, 1] = x_range[ix]
            x_data[ix, iy, iz, 2] = y_range[iy]
            x_data[ix, iy, iz, 3] = z_range[iz]
        end
        list = Batsrus.FileList("mock_3d.out", Batsrus.Real4Bat, ".", 0, 1, 0)
        bd = BATS(head, list, x_data, w_data)

        jx, jy, jz = get_current_density(bd)
        @test all(jx .== 0.0f0)
        @test all(jy .== 0.0f0)
        @test all(jz .== 0.0f0)
    end

    @testset "PLANETARY Units Scaling" begin
        nx, ny = 5, 5
        x_range = range(0.0f0, 4.0f0, length = nx)
        y_range = range(0.0f0, 4.0f0, length = ny)
        w_data = zeros(Float32, nx, ny, 3)
        # by = x
        for iy in 1:ny, ix in 1:nx
            w_data[ix, iy, 2] = x_range[ix]
        end

        head = Batsrus.BatsHead(
            2, SubString("PLANETARY"), 0, 0.0, false, 0, 3, [nx, ny], Float32[],
            [SubString("x"), SubString("y")],
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        # Provide real grid data to avoid division by zero
        x_data = zeros(Float32, nx, ny, 2)
        for iy in 1:ny, ix in 1:nx
            x_data[ix, iy, 1] = x_range[ix]
            x_data[ix, iy, 2] = y_range[iy]
        end
        bd = BATS(head, Batsrus.FileList("mock_p.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

        jx, jy, jz = get_current_density(bd)
        @test jz[2, 2] ≈ Float32(Batsrus.FAC_J_PLANETARY)
    end

    @testset "Raw Variable Fallback" begin
        nx, ny = 2, 2
        names = [
            SubString("bx"), SubString("by"), SubString("bz"),
            SubString("jx"), SubString("jy"), SubString("jz"),
        ]
        w_data = zeros(Float32, nx, ny, 6)
        w_data[:, :, 4] .= 1.23f0 # jx
        w_data[:, :, 5] .= 4.56f0 # jy
        w_data[:, :, 6] .= 7.89f0 # jz

        head = Batsrus.BatsHead(
            2, SubString("Mock"), 0, 0.0, false, 0, 6, [nx, ny], Float32[],
            [SubString("x"), SubString("y")], names, String[]
        )
        # Proper x_data even if not used directly for raw fallback
        x_data = zeros(Float32, nx, ny, 2)
        bd = BATS(head, Batsrus.FileList("mock_raw.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)

        @test bd[:jx][1, 1] == 1.23f0
        @test bd[:jy][1, 1] == 4.56f0
        @test bd[:jz][1, 1] == 7.89f0

        jx, jy, jz = get_current_density(bd)
        @test all(jx .== 1.23f0)
        @test all(jy .== 4.56f0)
        @test all(jz .== 7.89f0)
    end
end

end # module test_current_density
