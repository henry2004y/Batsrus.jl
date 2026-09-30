# Regression tests for the VTK AMR tree machinery that do not require data
# artifacts. The golden checksums were generated from the implementation that
# was verified equivalent to the original BATSRUS-derived code over an
# exhaustive case matrix (block sizes, all 27 neighbor positions, DiLevel
# 0-3, all siblings and root blocks) and against the 3D AMR connectivity
# SHA256 in tests_io.jl.

using Batsrus, StaticArrays, SHA
import Batsrus: Head, fillCellNeighbors!, parent_, child0_, unset_

# Build a fake Batl whose tree makes getSibling return `isibling` for node 5.
function make_batl(nI, nJ, nK, isibling, parent_set)
    iTree = fill(Int32(unset_), 18, 35)
    iTree[parent_, 5] = parent_set ? Int32(1) : Int32(unset_)
    # Give node 1 a full set of children; child `isibling` is node 5.
    for c in 1:8
        iTree[child0_ + c, 1] = c == isibling ? Int32(5) : Int32(10 + c)
    end
    head = Head(
        Int32(nI), Int32(nJ), Int32(nK), Int32(1), Int32(2), Int32(2), Int32(2),
        Int32(3), Int32(8), MVector{3, Int32}(1, 1, 1),
        MVector{3, Float64}(0, 0, 0), MVector{3, Float64}(1, 1, 1),
        Int32[1, 2, 3], MVector{3, Bool}(false, false, false),
        MVector{3, Float64}(0.1, 0.1, 0.1)
    )
    return Batl(head, iTree, MVector{3, Int32}(2, 2, 2), Int8(3))
end

@testset "fillCellNeighbors! golden" begin
    n2g = Int32.(1:35)
    ctx = SHA.SHA1_CTX()

    for (nI, nJ, nK) in ((4, 4, 4), (4, 4, 1))
        for di in 1:3, dj in 1:3, dk in 1:3   # all 27 neighbor positions
            for dlv in (Int8(-1), Int8(0), Int8(1), Int8(2), Int8(3))
                for s in 0:9, parent_set in (s == 0 ? (false,) : (true, false))
                    DiLevel = fill(Int8(-1), 3, 3, 3)
                    DiLevel[di, dj, dk] = dlv
                    iNode = fill(Int32(7), 4, 4, 4)
                    iNode[2, 2, 2] = Int32(5)      # this block's own node

                    batl = make_batl(nI, nJ, nK, max(s, 1), parent_set)
                    g = fill(Int32(-777), nI + 2, nJ + 2, nK + 2)
                    fillCellNeighbors!(batl, g, DiLevel, iNode, n2g)
                    SHA.update!(ctx, reinterpret(UInt8, [vec(g); nI; nJ; nK; di; dj; dk; dlv; s]))
                end
            end
        end
    end

    digest = bytes2hex(SHA.digest!(ctx))
    @test digest == "2c986b8e6cc8d8e24119363db0c5b942d694286e"
end
