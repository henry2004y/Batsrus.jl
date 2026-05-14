module BatsrusWriteVTKExt

using Batsrus, WriteVTK, Glob, LightXML
using PrecompileTools: @setup_workload, @compile_workload
using Batsrus: readtecdata, load, Batl, getConnectivity, get_range

"""
    convertTECtoVTU(file::AbstractString, outname="out")
    convertTECtoVTU(head, data, connectivity, outname="out")

Convert unstructured Tecplot data to VTK. Note that if using voxel type data in VTK, the
connectivity sequence is different from Tecplot: the 3D connectivity sequence in Tecplot is
the same as the `hexahedron` type in VTK, but different with the `voxel` type.
The 2D connectivity sequence is the same as the `quad` type, but different with the `pixel`
type. For example, in 3D the index conversion is:

```
# PLT to VTK voxel index_ = [1 2 4 3 5 6 8 7]
for i in 1:2
	connectivity = swaprows!(connectivity, 4*i-1, 4*i)
end
```
"""
function Batsrus.convertTECtoVTU(file::AbstractString, outname = "out")
    head, data, connectivity = readtecdata(file)
    return convertTECtoVTU(head, data, connectivity, outname)
end

function Batsrus.convertTECtoVTU(head, data, connectivity, outname = "out")
    nVar = length(head.variable)
    points = @view data[1:head.nDim, :]
    cells = Vector{MeshCell{VTKCellType, Array{Int32, 1}}}(undef, head.nCell)

    if head.nDim == 3
        @inbounds for i in 1:head.nCell
            cells[i] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, connectivity[:, i])
        end
    elseif head.nDim == 2
        @inbounds for i in 1:head.nCell
            cells[i] = MeshCell(VTKCellTypes.VTK_QUAD, connectivity[:, i])
        end
    end

    vtkfile = vtk_grid(outname, points, cells)

    for ivar in (head.nDim + 1):nVar
        if endswith(head.variable[ivar], "_x") # vector
            if head.nDim == 3
                var1 = @view data[ivar, :]
                var2 = @view data[ivar + 1, :]
                var3 = @view data[ivar + 2, :]
                namevar = replace(head.variable[ivar], "_x" => "")
                vtk_point_data(vtkfile, (var1, var2, var3), namevar)
            elseif head.nDim == 2
                var1 = @view data[ivar, :]
                var2 = @view data[ivar + 1, :]
                namevar = replace(head.variable[ivar], "_x" => "")
                vtk_point_data(vtkfile, (var1, var2), namevar)
            end
        elseif endswith(head.variable[ivar], r"_y|_z")
            continue
        else
            var = @view data[ivar, :]
            vtk_point_data(vtkfile, var, head.variable[ivar])
        end
    end

    # Add meta data from Tecplot AUXDATA
    for i in eachindex(head.auxdata)
        vtkfile[head.auxdataname[i], VTKFieldData()] = head.auxdata[i]
    end

    return outfiles = vtk_save(vtkfile)
end

"""
    convertIDLtoVTK(filename; gridType=1, verbose=false)

Convert 3D BATSRUS *.out to VTK. If `filename` does not end with "out", it tries to find the ".info" and ".tree" file with the same name tag and generates 3D unstructured VTU file.

# Keywords

  - `gridType::Symbol`: Type of target VTK grid (:vti: image, :vtr: rectilinear,
    :vts: structured grid, :vtm: multiblock, :vthb: AMR).
"""
function Batsrus.convertIDLtoVTK(
        filename::AbstractString;
        gridType::Symbol = :vti,
        outname = "out",
        verbose::Bool = false
    )
    if endswith(filename, ".out")
        bd = load(filename)
        nVar = length(bd.head.wname)

        if bd.head.ndim == 2
            if !bd.head.gencoord
                x, y = get_range(bd)
                outfiles = vtk_grid(outname, x, y) do vtk
                    for ivar in 1:nVar
                        if endswith(bd.head.wname[ivar], "x") # vector
                            var1 = @view bd.w[:, :, ivar]
                            var2 = @view bd.w[:, :, ivar + 1]
                            var3 = @view bd.w[:, :, ivar + 2]
                            namevar = bd.head.wname[ivar][1:(end - 1)]
                            vtk_point_data(vtk, (var1, var2, var3), namevar)
                        elseif endswith(bd.head.wname[ivar], "y") ||
                                endswith(bd.head.wname[ivar], "z")
                            continue
                        else
                            var = @view bd.w[:, :, ivar]
                            vtk_point_data(vtk, var, bd.head.wname[ivar])
                        end
                    end
                end
            else
                @error "2D generalized coordinate output conversion not supported!"
            end
        elseif bd.head.ndim == 3
            if gridType in (:vti, :vtr)
                if gridType == :vti # image
                    x, y, z = get_range(bd)
                else # rectilinear
                    x = @view bd.x[:, 1, 1, 1]
                    y = @view bd.x[1, :, 1, 2]
                    z = @view bd.x[1, 1, :, 3]
                end

                outfiles = vtk_grid(outname, x, y, z) do vtk
                    for ivar in 1:nVar
                        if bd.head.wname[ivar][end] == 'x' # vector
                            var1 = @view bd.w[:, :, :, ivar]
                            var2 = @view bd.w[:, :, :, ivar + 1]
                            var3 = @view bd.w[:, :, :, ivar + 2]
                            namevar = bd.head.wname[ivar][1:(end - 1)]
                            vtk_point_data(vtk, (var1, var2, var3), namevar)
                        elseif bd.head.wname[ivar][end] in ('y', 'z')
                            continue
                        else
                            var = @view bd.w[:, :, :, ivar]
                            vtk_point_data(vtk, var, bd.head.wname[ivar])
                        end
                    end
                end
            elseif gridType == :vts # structured grid
                xyz = permutedims(bd.x, [4, 1, 2, 3])

                outfiles = vtk_grid(outname, xyz) do vtk
                    for ivar in 1:nVar
                        if bd.head.wname[ivar][end] == 'x' # vector
                            var1 = @view bd.w[:, :, :, ivar]
                            var2 = @view bd.w[:, :, :, ivar + 1]
                            var3 = @view bd.w[:, :, :, ivar + 2]
                            namevar = bd.head.wname[ivar][1:(end - 1)]
                            vtk_point_data(vtk, (var1, var2, var3), namevar)
                        elseif bd.head.wname[ivar][end] in ('y', 'z')
                            continue
                        else
                            var = @view bd.w[:, :, :, ivar]
                            vtk_point_data(vtk, var, bd.head.wname[ivar])
                        end
                    end
                end
            else
                @error "No tree information for conversion!"
            end
        end
    else
        # info, tree, and out files
        bd = load(filename * ".out")
        batl = Batl(readhead(filename * ".info"), readtree(filename)...)

        if outname == "out"
            outname = filename
        end

        nDim = batl.nDim
        nVar = length(bd.head.wname)
        nI, nJ, nK = batl.head.nI, batl.head.nJ, batl.head.nK
        blockSize = nI * nJ * nK

        leaf_indices = findall(x -> x == Batsrus.used_,
            @view batl.iTree_IA[Batsrus.status_, :])
        # Sort leaf indices by (processor, local block index) to match .out file ordering
        sorted_leaf_indices = sort(leaf_indices,
            by = i -> (batl.iTree_IA[Batsrus.proc_, i], batl.iTree_IA[Batsrus.block_, i]))

        if gridType == :vthb
            # Native AMR support via vtkOverlappingAMR
            levels = @view batl.iTree_IA[Batsrus.level_, leaf_indices]
            unique_levels = sort(unique(levels))

            xdoc = XMLDocument()
            root = create_root(xdoc, "VTKFile")
            set_attributes(root; type = "vtkOverlappingAMR", version = "1.1",
                byte_order = "LittleEndian")

            amr = new_child(root, "vtkOverlappingAMR")
            set_attributes(amr;
                origin = "$(batl.head.CoordMin_D[1]) $(batl.head.CoordMin_D[2]) $(batl.head.CoordMin_D[3])",
                grid_description = "XYZ")

            # Spacing at level 0 (cell size)
            dx0 = (batl.head.CoordMax_D .- batl.head.CoordMin_D) ./
                  (batl.head.nRoot_D .* [nI, nJ, nK])

            outfiles = String[]

            for L in unique_levels
                block_node = new_child(amr, "Block")
                spacing_L = dx0 ./ (2^L)
                set_attributes(block_node; level = string(L),
                    spacing = "$(spacing_L[1]) $(spacing_L[2]) $(spacing_L[3])")

                level_leaf_indices = leaf_indices[levels .== L]

                for (idx_in_level, node_idx) in enumerate(level_leaf_indices)
                    # Find correct block index in the .out file
                    ib = findfirst(x -> x == node_idx, sorted_leaf_indices)
                    cell_range = (ib - 1) * blockSize + 1:ib * blockSize

                    c1 = batl.iTree_IA[Batsrus.coord1_, node_idx] - 1
                    c2 = batl.iTree_IA[Batsrus.coord2_, node_idx] - 1
                    c3 = size(batl.iTree_IA, 1) >= Batsrus.coord3_ ?
                         batl.iTree_IA[Batsrus.coord3_, node_idx] : Int32(1)
                    c3 = max(0, c3 - 1)
                    coord_D = [c1, c2, c3]

                    # Physical origin of this block
                    block_width_L = (batl.head.CoordMax_D .- batl.head.CoordMin_D) ./
                                    (batl.head.nRoot_D .* 2^L)
                    origin_block = batl.head.CoordMin_D .+ coord_D .* block_width_L

                    vti_name = "$(outname)_L$(L)_B$(idx_in_level)"

                    vtk = if nDim == 3
                        vtk_grid(vti_name, nI, nJ, nK; origin = Tuple(origin_block),
                            spacing = Tuple(spacing_L))
                    else
                        vtk_grid(vti_name, nI, nJ; origin = Tuple(origin_block[1:2]),
                            spacing = Tuple(spacing_L[1:2]))
                    end

                    for ivar in 1:nVar
                        if endswith(bd.head.wname[ivar], "x") # vector
                            if nDim == 3
                                var1 = reshape(@view(bd.w[cell_range, 1, 1, ivar]), nI, nJ, nK)
                                var2 = reshape(@view(bd.w[cell_range, 1, 1, ivar + 1]),
                                    nI, nJ, nK)
                                var3 = reshape(@view(bd.w[cell_range, 1, 1, ivar + 2]),
                                    nI, nJ, nK)
                                vtk_point_data(vtk, (var1, var2, var3),
                                    bd.head.wname[ivar][1:(end - 1)])
                            elseif nDim == 2
                                var1 = reshape(@view(bd.w[cell_range, 1, ivar]), nI, nJ)
                                var2 = reshape(@view(bd.w[cell_range, 1, ivar + 1]), nI, nJ)
                                vtk_point_data(vtk, (var1, var2),
                                    bd.head.wname[ivar][1:(end - 1)])
                            end
                        elseif endswith(bd.head.wname[ivar], r"y|z")
                            continue
                        else
                            if nDim == 3
                                var = reshape(@view(bd.w[cell_range, 1, 1, ivar]), nI, nJ, nK)
                            elseif nDim == 2
                                var = reshape(@view(bd.w[cell_range, 1, ivar]), nI, nJ)
                            end
                            vtk_point_data(vtk, var, bd.head.wname[ivar])
                        end
                    end
                    block_files = vtk_save(vtk)
                    append!(outfiles, block_files)

                    # Add to XML
                    ds = new_child(block_node, "DataSet")
                    lo = coord_D .* [nI, nJ, nK]
                    hi = (coord_D .+ 1) .* [nI, nJ, nK] .- 1
                    set_attributes(ds; index = string(idx_in_level - 1),
                        amr_box = "$(lo[1]) $(hi[1]) $(lo[2]) $(hi[2]) $(lo[3]) $(hi[3])",
                        file = splitdir(block_files[1])[2])
                end
            end
            vthb_name = outname * ".vthb"
            save_file(xdoc, vthb_name)
            pushfirst!(outfiles, vthb_name)
        elseif gridType == :vtm
            vtm = vtk_multiblock(outname)
            verbose && println("DEBUG: nLeaf = ", length(leaf_indices))

            for (i, node_idx) in enumerate(leaf_indices)
                # Find correct block index in the .out file
                ib = findfirst(x -> x == node_idx, sorted_leaf_indices)
                cell_range = (ib - 1) * blockSize + 1:ib * blockSize

                if nDim == 3
                    x_block = reshape(@view(bd.x[cell_range, 1, 1, 1]), nI, nJ, nK)
                    y_block = reshape(@view(bd.x[cell_range, 1, 1, 2]), nI, nJ, nK)
                    z_block = reshape(@view(bd.x[cell_range, 1, 1, 3]), nI, nJ, nK)

                    xi, yi, zi = x_block[:, 1, 1], y_block[1, :, 1], z_block[1, 1, :]

                    vtk_grid(vtm, xi, yi, zi) do vtk
                        for ivar in 1:nVar
                            if endswith(bd.head.wname[ivar], "x") # vector
                                var1 = reshape(@view(bd.w[cell_range, 1, 1, ivar]), nI, nJ, nK)
                                var2 = reshape(@view(bd.w[cell_range, 1, 1, ivar + 1]), nI, nJ, nK)
                                var3 = reshape(@view(bd.w[cell_range, 1, 1, ivar + 2]), nI, nJ, nK)
                                vtk_point_data(vtk, (var1, var2, var3), bd.head.wname[ivar][1:(end - 1)])
                            elseif endswith(bd.head.wname[ivar], r"y|z")
                                continue
                            else
                                var = reshape(@view(bd.w[cell_range, 1, 1, ivar]), nI, nJ, nK)
                                vtk_point_data(vtk, var, bd.head.wname[ivar])
                            end
                        end
                    end
                elseif nDim == 2
                    x_block = reshape(@view(bd.x[cell_range, 1, 1]), nI, nJ)
                    y_block = reshape(@view(bd.x[cell_range, 1, 2]), nI, nJ)

                    xi, yi = x_block[:, 1], y_block[1, :]

                    vtk_grid(vtm, xi, yi) do vtk
                        for ivar in 1:nVar
                            if endswith(bd.head.wname[ivar], "x") # vector
                                var1 = reshape(@view(bd.w[cell_range, 1, ivar]), nI, nJ)
                                var2 = reshape(@view(bd.w[cell_range, 1, ivar + 1]), nI, nJ)
                                vtk_point_data(vtk, (var1, var2), bd.head.wname[ivar][1:(end - 1)])
                            elseif endswith(bd.head.wname[ivar], r"y|z")
                                continue
                            else
                                var = reshape(@view(bd.w[cell_range, 1, ivar]), nI, nJ)
                                vtk_point_data(vtk, var, bd.head.wname[ivar])
                            end
                        end
                    end
                end
            end
            outfiles = vtk_save(vtm)
        else
            connectivity = getConnectivity(batl)
            nCell = size(connectivity, 2)
            if nDim == 3
                if batl.head.dxPlot_D[1] ≥ 0.0
                    @error "Why are there duplicate points?"
                else
                    points = bd.x[:, 1, 1, :]'
                end
            elseif nDim == 2
                if batl.head.dxPlot_D[1] ≥ 0.0 # points are not sorted in postproc.f90
                    @error "Why are there duplicate points? Ask!"
                    points = bd.x[:, 1, :]'
                else # points are sorted in postproc.f90!
                    @error "point original order cannot be retrieved!"
                end
            end
            cells = Vector{MeshCell{VTKCellType, Array{Int32, 1}}}(undef, nCell)

            if nDim == 3
                @inbounds for i in 1:nCell
                    cells[i] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, connectivity[:, i])
                end
            elseif nDim == 2
                @inbounds for i in 1:nCell
                    cells[i] = MeshCell(VTKCellTypes.VTK_QUAD, connectivity[:, i])
                end
            end

            vtkfile = vtk_grid(filename, points, cells)

            for ivar in 1:nVar
                if endswith(bd.head.wname[ivar], "x") # vector
                    if nDim == 3
                        var1 = @view bd.w[:, 1, 1, ivar]
                        var2 = @view bd.w[:, 1, 1, ivar + 1]
                        var3 = @view bd.w[:, 1, 1, ivar + 2]
                        var = (var1, var2, var3)
                        vtkfile[bd.head.wname[ivar][1:(end - 1)], VTKPointData()] = var
                    elseif nDim == 2 # not sure how VTK handles 2D vector!
                        var1 = @view bd.w[:, 1, ivar]
                        var2 = @view bd.w[:, 1, ivar + 1]
                        vtkfile[bd.head.wname[ivar][1:(end - 1)], VTKPointData()] = (var1, var2)
                    end
                elseif endswith(bd.head.wname[ivar], r"y|z")
                    continue
                else
                    if nDim == 3
                        var = @view bd.w[:, 1, 1, ivar]
                    elseif nDim == 2
                        var = @view bd.w[:, 1, ivar]
                    end
                    vtkfile[bd.head.wname[ivar], VTKPointData()] = var
                end
            end

            outfiles = vtk_save(vtkfile)
        end
    end

    verbose && @info "$(filename) finished conversion."

    return outfiles
end

"""
Return matrix X with swapped rows i and j.
"""
function swaprows!(X::Matrix, i::Int, j::Int)
    m, n = size(X)
    if (1 ≤ i ≤ m) && (1 ≤ j ≤ m)
        @inbounds @simd for k in 1:n
            X[i, k], X[j, k] = X[j, k], X[i, k]
        end
        return X
    else
        throw(BoundsError())
    end
end

"""
    create_pvd(filepattern)

Generate PVD file for a time series collection of VTK data.

# Example

```
create_pvd("*.vtu")
```
"""
function Batsrus.create_pvd(filepattern::String, debug::Bool = false)
    filenames = glob(filepattern)

    isempty(filenames) && @error "No files found matching pattern $(filepattern)!"

    # create an empty XML document
    xdoc = XMLDocument()

    # create & attach a root node
    xroot = create_root(xdoc, "VTKFile")

    type = "Collection"
    byte_order = "LittleEndian"
    compressor = "vtkZLibDataCompressor"

    set_attributes(xroot; type, byte_order, compressor)

    # create the first child
    xs1 = new_child(xroot, "Collection")

    for file in filenames
        i_end = findfirst("_n", file)[1] - 1

        second = parse(Int32, file[(i_end - 1):i_end])
        minute = parse(Int32, file[(i_end - 3):(i_end - 2)])
        timestep = 60 * minute + second
        debug && @show timestep

        e = new_child(xs1, "DataSet")

        set_attributes(e; timestep, group = "", part = "0", file)
    end

    name_index = findfirst("_t", filenames[1])[1] - 1

    save_file(xdoc, filenames[1][1:name_index] * ".pvd")

    return
end

@setup_workload begin
    @compile_workload begin
        mktempdir() do tmpdir
            # Mock data for precompilation
            file = joinpath(@__DIR__, "../test/precompile.out")
            outname = joinpath(tmpdir, "out")
            convertIDLtoVTK(file; outname)
        end
    end
end

end
