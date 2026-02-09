# All the IO related APIs.

# Fortran binary format includes record end/start tags that needs to be skipped.
# If there are continuous blocks, we usually skip 2*TAG between actual reading.
const TAG = 4 # Fortran record tag size

"""
    load(filename; npict=1, verbose=false)

Read BATSRUS output files. Stores the `npict` snapshot from an ascii or binary data file into the arrays of coordinates `x` and data `w`.
"""
function load(file::AbstractString; npict::Int = 1, verbose::Bool = false)
    filelist, fileID, pictsize = getfiletype(file)

    verbose && @info "filename=$(filelist.name)\n" * "npict=$(filelist.npictinfiles)"

    if filelist.npictinfiles - npict < 0
        throw(ArgumentError("Select snapshot $npict out of range $(filelist.npictinfiles)!"))
    end
    seekstart(fileID) # Rewind to start
    # Jump to the npict-th snapshot
    skip(fileID, pictsize * (npict - 1))

    head = getfilehead(fileID, filelist)

    # Read data
    if filelist.type == AsciiBat
        T = Float64
    else
        skip(fileID, TAG) # skip record start tag
        T = filelist.type == Real4Bat ? Float32 : Float64
    end

    x, w = allocateBuffer(head, T)
    _load_body!(x, w, fileID, filelist)

    close(fileID)

    return BATS(head, filelist, x, w)
end

function _load_body!(x::AbstractArray{T}, w::AbstractArray{T}, fileID::IOStream, filelist) where {T}
    return if filelist.type == AsciiBat
        getascii!(x, w, fileID)
    else
        getbinary!(x, w, fileID)
    end
end

"""
Read information from log file.
"""
function readlogdata(file::AbstractString)
    return open(file) do f
        nLine = countlines(f) - 2
        seekstart(f)
        headline = readline(f)
        variable = split(readline(f))
        ndim = 1
        it = 0
        t = 0.0
        gencoord = false
        nx = 1
        nw = length(variable)

        data = zeros(nw, nLine)
        @inbounds for i in 1:nLine
            line = split(readline(f))
            data[:, i] = Parsers.parse.(Float64, line)
        end

        head = (
            ndim = ndim, headline = headline, it = it, time = t, gencoord = gencoord,
            nw = nw, nx = nx, variable = variable,
        )

        return head, data
    end
end

"""
    readtecdata(file; verbose=false)

Return header, data and connectivity from BATSRUS Tecplot outputs. Both 2D and 3D binary and ASCII formats are supported.

# Examples

```
file = "3d_ascii.dat"
head, data, connectivity = readtecdata(file)
```
"""
function readtecdata(file::AbstractString; verbose::Bool = false)
    return open(file) do f
        head, pt0 = read_tecplot_header(f)

        data = Array{Float32, 2}(undef, length(head.variable), head.nNode)

        if head.nDim == 3
            connectivity = Array{Int32, 2}(undef, 8, head.nCell)
        elseif head.nDim == 2
            connectivity = Array{Int32, 2}(undef, 4, head.nCell)
        end

        # Check file type
        isBinary = false
        try
            Parsers.parse.(Float32, split(readline(f)))
        catch
            isBinary = true
            verbose && @info "reading binary file"
        end

        seek(f, pt0)

        if isBinary
            read_tecplot_data_binary!(f, data, connectivity)
        else
            read_tecplot_data_ascii!(f, data, connectivity)
        end

        return head, data, connectivity
    end
end

function read_tecplot_header(f)
    nDim = 3
    nNode = Int32(0)
    nCell = Int32(0)
    ET = ""
    title = ""
    VARS = SubString{String}[]

    # Read Tecplot header
    ln = readline(f) |> strip
    if startswith(ln, "TITLE")
        title = match(r"\"(.*?)\"", split(ln, '=', keepempty = false)[2])[1]
    else
        @warn "No title provided."
    end
    ln = readline(f) |> strip
    if startswith(ln, "VARIABLES")
        # Read until another keyword appears
        varline = split(ln, '=')[2]

        ln = readline(f) |> strip
        while !startswith(ln, "ZONE")
            varline *= strip(ln)
            ln = readline(f) |> strip
        end
        VARS = split(varline, '\"')
        deleteat!(VARS, findall(x -> x in (" ", ", ") || isempty(x), VARS))
    else
        @warn "No variable names provided."
    end

    while !startswith(ln, "AUXDATA")
        if !startswith(ln, "ZONE") # ZONE allows multiple \n
            zoneline = split(ln, ", ", keepempty = false)
        else # if the ZONE line has nothing, this won't work!
            zoneline = split(ln[6:end], ", ", keepempty = false)
            replace(zoneline[1], '"' => "") # Remove the quotes in T
        end
        for zline in zoneline
            name, value = split(zline, '=', keepempty = false)
            name = uppercase(name)
            if name == "T" # ZONE title
                T = value
            elseif name in ("NODES", "N")
                nNode = Parsers.parse(Int32, value)
            elseif name in ("ELEMENTS", "E")
                nCell = Parsers.parse(Int32, value)
            elseif name in ("ET", "ZONETYPE")
                if uppercase(value) in ("BRICK", "FEBRICK")
                    nDim = 3
                elseif uppercase(value) in ("QUADRILATERAL", "FEQUADRILATERAL")
                    nDim = 2
                end
                ET = uppercase(value)
            end
        end
        ln = readline(f) |> strip
    end

    auxdataname = String[]
    auxdata = Union{Int32, String}[]
    pt0 = position(f)

    while startswith(ln, "AUXDATA")
        name, value = split(ln, '"', keepempty = false)
        name = string(name[9:(end - 1)])
        str = string(strip(value))
        if name in ("ITER", "NPROC")
            str = Parsers.parse(Int32, value)
        elseif name == "TIMESIM"
            sec = split(str, "=")
            str = string(strip(sec[2]))
        end
        push!(auxdataname, name)
        push!(auxdata, str)

        pt0 = position(f)
        ln = readline(f) |> strip
    end

    if startswith(ln, "DT")
        pt0 = position(f)
    end

    seek(f, pt0)

    head = (
        variable = VARS, nNode = nNode, nCell = nCell, nDim = nDim, ET = ET,
        title = title, auxdataname = auxdataname, auxdata = auxdata,
    )

    return head, pt0
end

function read_tecplot_data_binary!(f, data, connectivity)
    nNode = size(data, 2)
    nCell = size(connectivity, 2)

    @inbounds for i in 1:nNode
        read!(f, @view data[:, i])
    end
    return @inbounds for i in 1:nCell
        read!(f, @view connectivity[:, i])
    end
end

function read_tecplot_data_ascii!(f, data::AbstractArray{T}, connectivity) where {T}
    nNode = size(data, 2)
    nCell = size(connectivity, 2)

    @inbounds for i in 1:nNode
        line = readline(f)
        offset = 1
        len = sizeof(line)
        bytes = codeunits(line)
        for j in axes(data, 1)
            # Skip local whitespace
            while offset <= len
                b = bytes[offset]
                if b == 0x20 || b == 0x09 || b == 0x0a || b == 0x0d
                    offset += 1
                else
                    break
                end
            end

            res = Parsers.xparse(T, line; pos = offset, len = len, delim = ' ')
            if (res.code & Parsers.OK) != 0
                data[j, i] = res.val
                offset += res.tlen
            end
        end
    end
    return @inbounds for i in 1:nCell
        line = readline(f)
        offset = 1
        len = sizeof(line)
        bytes = codeunits(line)
        for j in axes(connectivity, 1)
            # Skip local whitespace
            while offset <= len
                b = bytes[offset]
                if b == 0x20 || b == 0x09 || b == 0x0a || b == 0x0d
                    offset += 1
                else
                    break
                end
            end

            res = Parsers.xparse(Int32, line; pos = offset, len = len, delim = ' ')
            if (res.code & Parsers.OK) != 0
                connectivity[j, i] = res.val
                offset += res.tlen
            end
        end
    end
end

"""
Obtain file type.
"""
function getfiletype(file::AbstractString)
    fileID = open(file, "r")
    bytes = filesize(file)

    # Check the appendix of file names
    if occursin(r"^.*\.(log)$", file)
        type = LogBat
        npictinfiles = 1
    elseif occursin(r"^.*\.(dat)$", file) # Tecplot ascii format
        type = TecBat
        npictinfiles = 1
    else
        # Obtain filetype based on the length info in the first 4 bytes (Gabor's trick)
        lenhead = read(fileID, Int32)

        if lenhead != 79 && lenhead != 500
            type = AsciiBat
        else
            # The length of the 2nd line decides between real4 & real8
            # since it contains the time; which is real*8 | real*4
            skip(fileID, lenhead + TAG)
            len = read(fileID, Int32)
            if len == 20
                type = Real4Bat
            elseif len == 24
                type = Real8Bat
            else
                throw(ArgumentError("Incorrect formatted file: $file"))
            end
        end
        # Obtain file size & number of snapshots
        seekstart(fileID)
        pictsize = getfilesize(fileID, lenhead, Val(type))::Int
        npictinfiles = bytes ÷ pictsize
    end

    filelist = FileList(basename(file), type, dirname(file), bytes, npictinfiles, lenhead)

    return filelist, fileID, pictsize
end

"""
    getfilehead(fileID::IoStream, filelist::FileList) -> NameTuple

Obtain the header information from BATSRUS output file of `type` linked to `fileID`.

# Input arguments

  - `fileID::IOStream`: file identifier.
  - `filelist::FileList`: file information.
"""
function getfilehead(fileID::IOStream, filelist::FileList)
    if filelist.type == AsciiBat
        return _read_head_ascii(fileID)
    else # Real4Bat / Real8Bat
        return _read_head_binary(fileID, filelist.lenhead)
    end
end

function _read_head_ascii(fileID::IOStream)
    headline = readline(fileID)

    line_iter = eachsplit(readline(fileID))
    next = iterate(line_iter)
    it = Parsers.parse(Int32, next[1])
    next = iterate(line_iter, next[2])
    t = Parsers.parse(Float32, next[1])
    next = iterate(line_iter, next[2])
    ndim = Parsers.parse(Int32, next[1])
    next = iterate(line_iter, next[2])
    neqpar = Parsers.parse(Int32, next[1])
    next = iterate(line_iter, next[2])
    nw = Parsers.parse(Int32, next[1])

    gencoord = ndim < 0
    ndim = abs(ndim)

    nx = [Parsers.parse(Int32, s) for s in eachsplit(readline(fileID))]

    if neqpar > 0
        eqpar = [Parsers.parse(Float32, s) for s in eachsplit(readline(fileID))]
    else
        eqpar = Float32[]
    end

    varname = readline(fileID)

    return _create_batshead(ndim, headline, it, t, gencoord, neqpar, nw, nx, eqpar, varname)
end

function _read_head_binary(fileID::IOStream, lenstr::Integer)
    skip(fileID, TAG)
    headline = rstrip(String(read(fileID, lenstr)))
    skip(fileID, 2 * TAG)

    it = read(fileID, Int32)
    t = read(fileID, Float32)
    ndim = read(fileID, Int32)
    gencoord = ndim < 0
    ndim = abs(ndim)
    neqpar = read(fileID, Int32)
    nw = read(fileID, Int32)

    skip(fileID, 2 * TAG)
    nx = Vector{Int32}(undef, ndim)
    read!(fileID, nx)
    skip(fileID, 2 * TAG)

    if neqpar > 0
        eqpar = Vector{Float32}(undef, neqpar)
        read!(fileID, eqpar)
        skip(fileID, 2 * TAG)
    else
        eqpar = Float32[]
    end

    varname = String(read(fileID, lenstr))
    skip(fileID, TAG)

    return _create_batshead(ndim, headline, it, t, gencoord, neqpar, nw, nx, eqpar, varname)
end

function _create_batshead(ndim, headline, it, t, gencoord, neqpar, nw, nx, eqpar, varname)
    # Obtain output variable names
    variable = lowercase.(split(varname))
    coord = @view variable[1:ndim]
    wname = @view variable[(ndim + 1):(ndim + nw)]
    param = @view variable[(ndim + nw + 1):end]

    return BatsHead(ndim, headline, it, t, gencoord, neqpar, nw, nx, eqpar, coord, wname, param)
end

function skipline(s::IO)
    while !eof(s)
        c = read(s, Char)
        c == '\n' && break
    end

    return
end

"""
    getfilesize(fileID::IOStream, lenstr::Int32, ::Val{FileType})

Return the size in bytes for one snapshot.
"""
function _getfilesize_binary(fileID::IOStream, lenstr::Integer, ::Val{T}) where {T}
    pointer0 = position(fileID) # Record header start location

    skip(fileID, TAG + lenstr + 2 * TAG + sizeof(Int32) + sizeof(Float32))
    ndim = abs(read(fileID, Int32))
    nt = read(fileID, Int32)
    nw = read(fileID, Int32)
    skip(fileID, 2 * TAG)

    prod_nx = 1
    for _ in 1:ndim
        prod_nx *= read(fileID, Int32)
    end

    skip(fileID, 2 * TAG)
    if nt > 0
        skip(fileID, nt * sizeof(Float32))
        skip(fileID, 2 * TAG)
    end
    skip(fileID, lenstr)
    skip(fileID, TAG)

    pointer1 = position(fileID)
    headlen = pointer1 - pointer0 # header length

    type_size = T === Real4Bat ? 4 : 8

    # Calculate the snapshot size = header + data + recordmarks
    # 8 bytes = 2 * TAG (4 bytes each)
    return pictsize = headlen + 8 * (1 + nw) + type_size * (ndim + nw) * prod_nx
end

function getfilesize(fileID::IOStream, lenstr::Int32, v::Val{Real4Bat})
    return _getfilesize_binary(fileID, lenstr, v)
end
function getfilesize(fileID::IOStream, lenstr::Int32, v::Val{Real8Bat})
    return _getfilesize_binary(fileID, lenstr, v)
end

function getfilesize(fileID::IOStream, lenstr::Int32, ::Val{AsciiBat})
    pointer0 = position(fileID) # Record header start location

    skipline(fileID)
    line = readline(fileID)
    line = split(line)
    ndim = Parsers.parse(Int32, line[3])
    neqpar = Parsers.parse(Int32, line[4])
    nw = Parsers.parse(Int32, line[5])
    ndim = abs(ndim)
    nx = Parsers.parse.(Int64, split(readline(fileID)))
    neqpar > 0 && skipline(fileID)
    skipline(fileID)

    pointer1 = position(fileID)
    headlen = pointer1 - pointer0 # header length
    # Calculate the snapshot size = header + data + recordmarks
    return pictsize = headlen + (18 * (ndim + nw) + 1) * prod(nx)
end

getfilesize(fileID::IOStream, lenstr::Int32, ::Val{LogBat}) = 1

"""
Create buffer for x and w.
"""
allocateBuffer(head::BatsHead, ::Type{T}) where {T} = _allocateBuffer(Val(head.ndim), head, T)

function _allocateBuffer(::Val{N}, head::BatsHead, ::Type{T}) where {N, T}
    dims = ntuple(i -> head.nx[i], Val(N))
    x = Array{T, N + 1}(undef, dims..., N)
    w = Array{T, N + 1}(undef, dims..., head.nw)
    return x, w
end

"""
Read ascii format coordinates and data values.
"""
function getascii!(x::AbstractArray{T, N}, w::AbstractArray{T, N}, fileID::IOStream) where {T, N}
    # N is total dimensions (including species/components).
    # The spatial dimensions are 1:N-1.
    spatial_dims = ntuple(i -> size(x, i), Val(N - 1))
    ndim = size(x, N)
    nw = size(w, N)

    return @inbounds for id in CartesianIndices(spatial_dims)
        line = readline(fileID)
        offset = 1
        len = sizeof(line)
        bytes = codeunits(line)

        # Read coordinates
        for k in 1:ndim
            # Skip local whitespace
            while offset <= len
                b = bytes[offset]
                if b == 0x20 || b == 0x09 || b == 0x0a || b == 0x0d
                    offset += 1
                else
                    break
                end
            end

            res = Parsers.xparse(T, line; pos = offset, len = len, delim = ' ')
            if (res.code & Parsers.OK) != 0
                x[id, k] = res.val
                offset += res.tlen
            end
        end

        # Read variables
        for k in 1:nw
            # Skip local whitespace
            while offset <= len
                b = bytes[offset]
                if b == 0x20 || b == 0x09 || b == 0x0a || b == 0x0d
                    offset += 1
                else
                    break
                end
            end

            res = Parsers.xparse(T, line; pos = offset, len = len, delim = ' ')
            if (res.code & Parsers.OK) != 0
                w[id, k] = res.val
                offset += res.tlen
            end
        end
    end
end

"""
Read binary format coordinates and data values.
"""
function getbinary!(x::AbstractArray{T, N}, w, fileID::IOStream) where {T, N}
    read!(fileID, x)
    skip(fileID, 2 * TAG)
    return @inbounds for iw in axes(w, N)
        read!(fileID, selectdim(w, N, iw))
        skip(fileID, 2 * TAG)
    end
end

function Base.show(io::IO, data::BatsrusIDL)
    showhead(io, data)
    if data.list.bytes ≥ 1.0e9
        str = @sprintf "filesize: %.1f GB" data.list.bytes / 1.0e9
    elseif data.list.bytes ≥ 1.0e6
        str = @sprintf "filesize: %.1f MB" data.list.bytes / 1.0e6
    elseif data.list.bytes ≥ 1.0e3
        str = @sprintf "filesize: %.1f KB" data.list.bytes / 1.0e3
    else
        str = @sprintf "filesize: %.1f Bytes" data.list.bytes
    end
    println(io, str)
    return println(io, "snapshots: ", data.list.npictinfiles)
end

"""
    showhead(file, head)

Displaying file header information.
"""
function showhead(file::FileList, head::BatsHead, io::IO = stdout)
    print(io, "filename : ")
    printstyled(io, file.name, '\n'; color = :cyan, underline = true)
    println(io, "filetype : ", file.type)
    println(io, "headline : ", head.headline)
    print(io, "iteration: ")
    printstyled(io, head.it, '\n'; color = :cyan)
    print(io, "time     : ")
    printstyled(io, head.time, '\n'; color = :cyan)
    print(io, "gencoord: ")
    printstyled(io, head.gencoord, '\n'; color = :yellow)
    print(io, "ndim     : ")
    printstyled(io, head.ndim, '\n'; color = :yellow)

    if head.neqpar > 0
        print(io, "parameters: [ ")
        for par in head.eqpar
            print(io, par, " ")
        end
        print(io, "]\ncoordinate names: [ ")
        for c in head.coord
            print(io, c, " ")
        end
        print(io, "]\nvariable   names: [ ")
        for w in head.wname
            printstyled(io, w, " "; color = :green)
        end
        print(io, "]\nparameter  names: [ ")
        for p in head.param
            print(io, p, " ")
        end
        println(io, "]")
    end

    return
end

"""
    showhead(data)
    showhead(io, data)

Display file information of `data`.
"""
showhead(data::BatsrusIDL) = showhead(data.list, data.head)
showhead(io::IO, data::BatsrusIDL) = showhead(data.list, data.head, io)
