module BatsrusWriteVTKExt

using Batsrus, WriteVTK, Glob, LightXML
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

  - `gridType::Symbol`: Type of target VTK grid (vti: image, vtr: rectilinear, vts: structured grid).
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
                        if bd.head.wname[ivar][end] == 'x' # vector
                            var1 = @view bd.w[:, :, ivar]
                            var2 = @view bd.w[:, :, ivar + 1]
                            var3 = @view bd.w[:, :, ivar + 2]
                            namevar = bd.head.wname[ivar][1:(end - 1)]
                            vtk_point_data(vtk, (var1, var2, var3), namevar)
                        elseif bd.head.wname[ivar][end] in ('y', 'z')
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
        connectivity = getConnectivity(batl)

        outname = filename

        nDim = batl.nDim
        nVar = length(bd.head.wname)
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

    verbose && @info "$(filename) finished conversion."

    return
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
function Batsrus.create_pvd(filepattern::String)
    filenames = glob(filepattern)

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
        @show timestep

        e = new_child(xs1, "DataSet")

        set_attributes(e; timestep, group = "", part = "0", file)
    end

    name_index = findfirst("_t", filenames[1])[1] - 1

    save_file(xdoc, filenames[1][1:name_index] * ".pvd")

    return
end

end
