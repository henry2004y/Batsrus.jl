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
        
        head = Batsrus.BatsHead(
            1, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx], Float32[],
            [SubString("x")], [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        x_data = collect(reshape(x_range, nx, 1))
        list = Batsrus.FileList("mock_1d.out", Batsrus.Real4Bat, ".", 0, 1, 0)
        bd = BATS(head, list, x_data, w_data)
        
        jx, jy, jz = get_current_density(bd)
        @test all(jx .== 0.0f0)
        @test all(jy .== -0.0f0)
        @test jz[2:4] ≈ [1.0f0, 1.0f0, 1.0f0]
        
        # Test individual component access (triggers _compute_j*)
        @test all(bd[:jx] .== 0.0f0)
        @test all(bd[:jy] .== 0.0f0)
        @test bd[:jz][2:4] ≈ [1.0f0, 1.0f0, 1.0f0]
    end

    @testset "2D Structured" begin
        # B_z = x, so J_y = -dBz/dx = -1
        # B_z = y, so J_x = dBz/dy = 1
        nx, ny = 5, 5
        x_range = range(0.0f0, 4.0f0, length = nx)
        y_range = range(0.0f0, 4.0f0, length = ny)
        
        w_data = zeros(Float32, nx, ny, 3)
        for iy in 1:ny, ix in 1:nx
            w_data[ix, iy, 3] = x_range[ix] + y_range[iy] # bz = x + y
        end
        
        head = Batsrus.BatsHead(
            2, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx, ny], Float32[],
            [SubString("x"), SubString("y")], 
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        # BATS expects (nx, ny, 2) for x_data in 2D
        x_data = zeros(Float32, nx, ny, 2)
        for iy in 1:ny, ix in 1:nx
            x_data[ix, iy, 1] = x_range[ix]
            x_data[ix, iy, 2] = y_range[iy]
        end
        bd = BATS(head, Batsrus.FileList("mock_2d.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)
        
        jx, jy, jz = get_current_density(bd)
        @test jx[2:4, 2:4] ≈ ones(Float32, 3, 3)    # dBz/dy = 1
        @test jy[2:4, 2:4] ≈ -ones(Float32, 3, 3)   # -dBz/dx = -1
        @test all(jz .== 0.0f0)

        # Test individual component access
        @test bd[:jx][2:4, 2:4] ≈ ones(Float32, 3, 3)
        @test bd[:jy][2:4, 2:4] ≈ -ones(Float32, 3, 3)
        @test all(bd[:jz] .== 0.0f0)
    end

    @testset "3D Structured" begin
        nx, ny, nz = 5, 5, 5
        x_range = range(0.0f0, 4.0f0, length = nx)
        y_range = range(0.0f0, 4.0f0, length = ny)
        z_range = range(0.0f0, 4.0f0, length = nz)
        
        w_data = zeros(Float32, nx, ny, nz, 3)
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            w_data[ix, iy, iz, 1] = y_range[iy] # bx = y
            w_data[ix, iy, iz, 2] = x_range[ix] # by = x
        end
        
        head = Batsrus.BatsHead(
            3, SubString("Mock"), 0, 0.0f0, false, 0, 3, [nx, ny, nz], Float32[],
            [SubString("x"), SubString("y"), SubString("z")], 
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        x_data = zeros(Float32, nx, ny, nz, 3)
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            x_data[ix, iy, iz, 1] = x_range[ix]
            x_data[ix, iy, iz, 2] = y_range[iy]
            x_data[ix, iy, iz, 3] = z_range[iz]
        end
        bd = BATS(head, Batsrus.FileList("mock_3d.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)
        
        jx, jy, jz = get_current_density(bd)
        @test all(jx .== 0.0f0)
        @test all(jy .== 0.0f0)
        @test all(isapprox.(jz[2:4, 2:4, 2:4], 0.0f0, atol=1e-6))
        
        # Test individual component access
        @test all(bd[:jx] .== 0.0f0)
        @test all(bd[:jy] .== 0.0f0)
        @test all(isapprox.(bd[:jz][2:4, 2:4, 2:4], 0.0f0, atol=1e-6))
    end

    @testset "Unstructured/Curvilinear Grid Error" begin
        # For ndim=1, BATS expects (nx, 1) for x_data
        nx = 4
        head = Batsrus.BatsHead(
            1, SubString("Mock"), 0, 0.0f0, true, 0, 3, [nx], Float32[],
            [SubString("x")], 
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        x_data = zeros(Float32, nx, 2) # Dimp1 = 2 = ndim+1
        w_data = zeros(Float32, nx, 3)
        bd = BATS(head, Batsrus.FileList("mock_u.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)
        
        @test_throws Exception get_current_density(bd)
    end

    @testset "PLANETARY Units Scaling" begin
        nx, ny = 5, 5
        x_range = range(0.0f0, 4.0f0, length = nx)
        y_range = range(0.0f0, 4.0f0, length = ny)
        w_data = zeros(Float32, nx, ny, 3)
        for iy in 1:ny, ix in 1:nx
            w_data[ix, iy, 2] = x_range[ix]
        end
        
        head = Batsrus.BatsHead(
            2, SubString("PLANETARY"), 0, 0.0f0, false, 0, 3, [nx, ny], Float32[],
            [SubString("x"), SubString("y")], 
            [SubString("bx"), SubString("by"), SubString("bz")], String[]
        )
        x_data = zeros(Float32, nx, ny, 2)
        for iy in 1:ny, ix in 1:nx
            x_data[ix, iy, 1] = x_range[ix]
            x_data[ix, iy, 2] = y_range[iy]
        end
        bd = BATS(head, Batsrus.FileList("mock_p.out", Batsrus.Real4Bat, ".", 0, 1, 0), x_data, w_data)
        
        jx, jy, jz = get_current_density(bd)
        @test jz[2, 2] ≈ Float32(Batsrus.FAC_J_PLANETARY)
    end
end

end # module test_current_density
