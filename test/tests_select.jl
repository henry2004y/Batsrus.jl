
using Batsrus
using DimensionalData
using StaticArrays
using Test

# Mock BATS creation
function create_mock_bats(gencoord::Bool)
   nx, ny, nz = 10, 10, 10
   ndim = 3
   nw = 1

   # Coordinates
   xrange = range(0.0, 1.0, length = nx)
   yrange = range(0.0, 1.0, length = ny)
   zrange = range(0.0, 1.0, length = nz)

   x_data = zeros(Float32, nx, ny, nz, 3)
   for (i, x) in enumerate(xrange), (j, y) in enumerate(yrange), (k, z) in enumerate(zrange)
      x_data[i, j, k, 1] = x
      x_data[i, j, k, 2] = y
      x_data[i, j, k, 3] = z
   end

   # Variable data
   w_data = zeros(Float32, nx, ny, nz, nw)
   for i in 1:nx, j in 1:ny, k in 1:nz
      w_data[i, j, k, 1] = i + j + k
   end

   head = Batsrus.BatsHead(
      ndim, "Mock Data $(gencoord)", 0, 0.0, gencoord, 0, nw, [nx, ny, nz], Float32[],
      ["x", "y", "z"], ["var1"], String[]
   )

   filetype = Batsrus.Real4Bat
   list = Batsrus.FileList(
      "mock.out", filetype, ".", 0, 1, 0
   )

   return BATS(head, list, x_data, w_data)
end

@testset "Select Functions" begin
   for gencoord in [false, true]
      @testset "gencoord=$gencoord" begin
         bd = create_mock_bats(gencoord)

         # Test cutdata
         # Basic cut without plotrange should works for both (using selectdim)
         cut1, cut2, W = cutdata(bd, "var1"; dir = "x", sequence = 5)
         @test size(W) == (10, 10)

         if gencoord
            # gencoord=true: Expect error if we try to use subsurface/subvolume via plotrange
            @test_throws ErrorException cutdata(
               bd, "var1"; dir = "x", sequence = 5, plotrange = [0.2, 0.8, 0.2, 0.8])

            # Verify explicit subvolume call fails
            limits_vol = [0.2, 0.8, 0.2, 0.8, 0.2, 0.8]
            @test_throws ErrorException subvolume(
               bd.x[:, :, :, 1], bd.x[:, :, :, 2], bd.x[:, :, :, 3],
               bd.w[:, :, :, 1], limits_vol
            )
         else
            # gencoord=false
            # Test cutdata with plotrange (calls subsurface)
            cut1_sub, cut2_sub, W_sub = cutdata(
               bd, "var1"; dir = "x", sequence = 5, plotrange = [0.2, 0.8, 0.2, 0.8])

            # Check dimensions. Range 0.0 to 1.0 with 10 points:
            # 0.0, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 1.0
            # [0.2, 0.8] should include: 0.22, 0.33, 0.44, 0.55, 0.66, 0.77 (6 points)
            @test size(W_sub) == (6, 6)

            # Test explicit subsurface
            limits = [0.2, 0.8, 0.2, 0.8]
            subx, suby, subW = subsurface(cut1, cut2, W, limits)
            @test size(subW) == (6, 6)

            # Test subvolume
            limits_vol = [0.2, 0.8, 0.2, 0.8, 0.2, 0.8]
            subx_v, suby_v, subz_v, subW_v = subvolume(
               bd.x[:, :, :, 1], bd.x[:, :, :, 2], bd.x[:, :, :, 3],
               bd.w[:, :, :, 1], limits_vol
            )
            @test size(subW_v) == (6, 6, 6)
         end
      end
   end
end
