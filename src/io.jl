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

   verbose && @info "filename=$(filelist.name)\n"*"npict=$(filelist.npictinfiles)"

   if filelist.npictinfiles - npict < 0
      throw(ArgumentError("Select snapshot $npict out of range $(filelist.npictinfiles)!"))
   end
   seekstart(fileID) # Rewind to start
   # Jump to the npict-th snapshot
   skip(fileID, pictsize*(npict-1))

   head = getfilehead(fileID, filelist)
   # Read data
   if filelist.type == AsciiBat
      T = Float64 # why Float64?
      x, w = allocateBuffer(head, T)
      getascii!(x, w, fileID)
   else
      skip(fileID, TAG) # skip record start tag
      T = filelist.type == Real4Bat ? Float32 : Float64
      x, w = allocateBuffer(head, T)
      getbinary!(x, w, fileID)
   end

   close(fileID)

   #setunits(head,"")

   BATS(head, filelist, x, w)
end

"""
Read information from log file.
"""
function readlogdata(file::AbstractString)
   f = open(file)
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

   close(f)

   head = (ndim = ndim, headline = headline, it = it, time = t, gencoord = gencoord,
      nw = nw, nx = nx, variable = variable)

   head, data
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
   f = open(file)

   nDim = 3
   nNode = Int32(0)
   nCell = Int32(0)
   ET = ""

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
         replace(zoneline[1], '"'=>"") # Remove the quotes in T
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

   while startswith(ln, "AUXDATA") || startswith(ln, "DT")
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

   seek(f, pt0)

   data = Array{Float32, 2}(undef, length(VARS), nNode)

   if nDim == 3
      connectivity = Array{Int32, 2}(undef, 8, nCell)
   elseif nDim == 2
      connectivity = Array{Int32, 2}(undef, 4, nCell)
   end

   isBinary = false
   try
      Parsers.parse.(Float32, split(readline(f)))
   catch
      isBinary = true
      verbose && @info "reading binary file"
   end

   seek(f, pt0)

   if isBinary
      @inbounds for i in 1:nNode
         read!(f, @view data[:, i])
      end
      @inbounds for i in 1:nCell
         read!(f, @view connectivity[:, i])
      end
   else
      @inbounds for i in 1:nNode
         x = readline(f)
         data[:, i] .= Parsers.parse.(Float32, split(x))
      end
      @inbounds for i in 1:nCell
         x = readline(f)
         connectivity[:, i] .= Parsers.parse.(Int32, split(x))
      end
   end

   close(f)

   head = (variable = VARS, nNode = nNode, nCell = nCell, nDim = nDim, ET = ET,
      title = title, auxdataname = auxdataname, auxdata = auxdata)

   head, data, connectivity
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

   filelist, fileID, pictsize
end

"""
     getfilehead(fileID::IoStream, filelist::FileList) -> NameTuple

Obtain the header information from BATSRUS output file of `type` linked to `fileID`.

# Input arguments

  - `fileID::IOStream`: file identifier.
  - `filelist::FileList`: file information.
"""
function getfilehead(fileID::IOStream, filelist::FileList)
   type, lenstr = filelist.type, filelist.lenhead

   if type == AsciiBat
      headline = readline(fileID)
      line = readline(fileID) |> split
      it = Parsers.parse(Int32, line[1])
      t = Parsers.parse(Float32, line[2])
      ndim = Parsers.parse(Int32, line[3])
      neqpar = Parsers.parse(Int32, line[4])
      nw = Parsers.parse(Int32, line[5])
      gencoord = ndim < 0
      ndim = abs(ndim)
      nx = Parsers.parse.(Int32, split(readline(fileID)))
      if neqpar > 0
         eqpar = Parsers.parse.(Float32, split(readline(fileID)))
      end
      varname = readline(fileID)
   elseif type ∈ (Real4Bat, Real8Bat)
      skip(fileID, TAG)
      headline = rstrip(String(read(fileID, lenstr)))
      skip(fileID, 2*TAG)
      it = read(fileID, Int32)
      t = read(fileID, Float32)
      ndim = read(fileID, Int32)
      gencoord = ndim < 0
      ndim = abs(ndim)
      neqpar = read(fileID, Int32)
      nw = read(fileID, Int32)
      skip(fileID, 2*TAG)
      nx = zeros(Int32, ndim)
      read!(fileID, nx)
      skip(fileID, 2*TAG)
      if neqpar > 0
         eqpar = zeros(Float32, neqpar)
         read!(fileID, eqpar)
         skip(fileID, 2*TAG)
      end
      varname = String(read(fileID, lenstr))
      skip(fileID, TAG)
   end

   # Obtain output variable names
   variable = lowercase.(split(varname))
   coord = @view variable[1:ndim]
   wname = @view variable[(ndim + 1):(ndim + nw)]
   param = @view variable[(ndim + nw + 1):end]

   head = BatsHead(ndim, headline, it, t, gencoord, neqpar, nw, nx, eqpar,
      coord, wname, param)
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
function getfilesize(fileID::IOStream, lenstr::Int32, ::Val{Real4Bat})
   pointer0 = position(fileID) # Record header start location

   skip(fileID, TAG + lenstr + 2*TAG + sizeof(Int32) + sizeof(Float32))
   ndim = abs(read(fileID, Int32))
   nt = read(fileID, Int32)
   nw = read(fileID, Int32)
   skip(fileID, 2*TAG)
   nx = Vector{Int32}(undef, ndim)
   read!(fileID, nx)
   skip(fileID, 2*TAG)
   if nt > 0
      tmp = zeros(Float32, nt)
      read!(fileID, tmp)
      skip(fileID, 2*TAG)
   end
   read(fileID, lenstr)
   skip(fileID, TAG)

   pointer1 = position(fileID)
   headlen = pointer1 - pointer0 # header length
   # Calculate the snapshot size = header + data + recordmarks
   pictsize = headlen + 8*(1 + nw) + 4*(ndim + nw)*prod(nx)
end

function getfilesize(fileID::IOStream, lenstr::Int32, ::Val{Real8Bat})
   pointer0 = position(fileID) # Record header start location

   skip(fileID, TAG + lenstr + 2*TAG + sizeof(Int32) + sizeof(Float32))
   ndim = abs(read(fileID, Int32))
   nt = read(fileID, Int32)
   nw = read(fileID, Int32)
   skip(fileID, 2*TAG)
   nx = Vector{Int32}(undef, ndim)
   read!(fileID, nx)
   skip(fileID, 2*TAG)
   if nt > 0
      tmp = zeros(Float32, nt)
      read!(fileID, tmp)
      skip(fileID, 2*TAG)
   end
   read(fileID, lenstr)
   skip(fileID, TAG)

   pointer1 = position(fileID)
   headlen = pointer1 - pointer0 # header length
   # Calculate the snapshot size = header + data + recordmarks
   pictsize = headlen + 8*(1 + nw) + 8*(ndim + nw)*prod(nx)
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
   pictsize = headlen + (18*(ndim + nw) + 1)*prod(nx)
end

getfilesize(fileID::IOStream, lenstr::Int32, ::Val{LogBat}) = 1

"""
Create buffer for x and w.
"""
function allocateBuffer(head::BatsHead, T::DataType)
   if head.ndim == 1
      n1 = head.nx[1]
      x = Array{T, 2}(undef, n1, head.ndim)
      w = Array{T, 2}(undef, n1, head.nw)
   elseif head.ndim == 2
      n1, n2 = head.nx
      x = Array{T, 3}(undef, n1, n2, head.ndim)
      w = Array{T, 3}(undef, n1, n2, head.nw)
   elseif head.ndim == 3
      n1, n2, n3 = head.nx
      x = Array{T, 4}(undef, n1, n2, n3, head.ndim)
      w = Array{T, 4}(undef, n1, n2, n3, head.nw)
   end

   x, w
end

"""
Read ascii format coordinates and data values.
"""
function getascii!(x::Array{T, 2}, w, fileID::IOStream) where T
   @inbounds @views for id in axes(x, 1)
      temp = Parsers.parse.(Float64, split(readline(fileID)))
      x[id] = temp[1]
      w[id, :] = temp[2:end]
   end
end

function getascii!(x::Array{T, 3}, w, fileID::IOStream) where T
   @inbounds @views for id in CartesianIndices(size(x)[1:2])
      temp = Parsers.parse.(Float64, split(readline(fileID)))
      x[id, :] = temp[1:2]
      w[id, :] = temp[3:end]
   end
end

function getascii!(x::Array{T, 4}, w, fileID::IOStream) where T
   @inbounds @views for id in CartesianIndices(size(x)[1:3])
      temp = Parsers.parse.(Float64, split(readline(fileID)))
      x[id, :] = temp[1:3]
      w[id, :] = temp[4:end]
   end
end

"""
Read binary format coordinates and data values.
"""
function getbinary!(x::Array{T, 2}, w, fileID::IOStream) where T
   read!(fileID, x)
   skip(fileID, 2*TAG)
   dimlast = 2
   @inbounds for iw in axes(w, dimlast)
      read!(fileID, selectdim(w, dimlast, iw))
      skip(fileID, 2*TAG)
   end
end

function getbinary!(x::Array{T, 3}, w, fileID::IOStream) where T
   read!(fileID, x)
   skip(fileID, 2*TAG)
   dimlast = 3
   @inbounds for iw in axes(w, dimlast)
      read!(fileID, selectdim(w, dimlast, iw))
      skip(fileID, 2*TAG)
   end
end

function getbinary!(x::Array{T, 4}, w, fileID::IOStream) where T
   read!(fileID, x)
   skip(fileID, 2*TAG)
   dimlast = 4
   @inbounds for iw in axes(w, dimlast)
      read!(fileID, selectdim(w, dimlast, iw))
      skip(fileID, 2*TAG)
   end
end

"""
     setunits(head, type; distance=1.0, mp=1.0, me=1.0)

Set the units for the output files.
If type is given as "SI", "CGS", "NORMALIZED", "PIC", "PLANETARY", "SOLAR", set
`typeunit = type`, otherwise try to guess from the fileheader. Based on `typeunit` set units
for distance [xSI], time [tSI], density [rhoSI], pressure [pSI], magnetic field [bSI] and
current density [jSI] in SI units. Distance unit [rplanet | rstar], ion and electron mass in
[amu] can be set with optional `distance`, `mp` and `me`.

Also calculate convenient constants ti0, cs0 ... for typical formulas.
This function is currently not used anywhere!
"""
function setunits(head::BatsHead, type; distance = 1.0, mp = 1.0, me = 1.0)
   ndim = head.ndim
   headline = head.headline
   neqpar = head.neqpar
   nw = head.nw
   eqpar = head.eqpar
   param = head.param

   mu0SI = 4π*1e-7      # H/m
   cSI = 2.9978e8       # speed of light, [m/s]
   mpSI = 1.6726e-27     # kg
   eSI = 1.602e-19      # elementary charge, [C]
   AuSI = 149597870700   # m
   RsSI = 6.957e8        # m
   Mi = 1.0            # Ion mass, [amu]
   Me = 1.0/1836.15    # Electron mass, [amu]
   gamma = 5/3            # Adiabatic index for first fluid
   gammae = 5/3            # Adiabatic index for electrons
   kbSI = 1.38064852e-23 # Boltzmann constant, [m2 kg s-2 K-1]
   e0SI = 8.8542e-12     # [F/m]

   # This part is used to guess the units.
   # To be honest; I don`t understand the logic here. For example
   # nPa & m/s may appear in the same headline?
   if type !== ""
      typeunit = uppercase(type)
   elseif occursin("PIC", headline)
      typeunit = "PIC"
   elseif occursin(" AU ", headline)
      typeunit = "OUTERHELIO"
   elseif occursin(r"(kg/m3)|(m/s)", headline)
      typeunit = "SI"
   elseif occursin(r"(nPa)|( nT )", headline)
      typeunit = "PLANETARY"
   elseif occursin(r"(dyne)|( G)", headline)
      typeunit = "SOLAR"
   else
      typeunit = "NORMALIZED"
   end

   if typeunit == "SI"
      xSI = 1.0             # m
      tSI = 1.0             # s
      rhoSI = 1.0             # kg/m^3
      uSI = 1.0             # m/s
      pSI = 1.0             # Pa
      bSI = 1.0             # T
      jSI = 1.0             # A/m^2
   elseif typeunit == "CGS"
      xSI = 0.01            # cm
      tSI = 1.0             # s
      rhoSI = 1000.0          # g/cm^3
      uSI = 0.01            # cm/s
      pSI = 0.1             # dyne/cm^2
      bSI = 1.0e-4          # G
      jSI = 10*cSI          # Fr/s/cm^2
   elseif typeunit == "PIC"
      # Normalized PIC units
      xSI = 1.0             # cm
      tSI = 1.0             # s
      rhoSI = 1.0             # g/cm^3
      uSI = 1.0             # cm/s
      pSI = 1.0             # dyne/cm^2
      bSI = 1.0             # G
      jSI = 1.0             # Fr/s/cm^2
      c0 = 1.0             # speed of light always 1 for iPIC3D
   elseif typeunit == "NORMALIZED"
      xSI = 1.0             # distance unit in SI
      tSI = 1.0             # time unit in SI
      rhoSI = 1.0             # density unit in SI
      uSI = 1.0             # velocity unit in SI
      pSI = 1.0             # pressure unit in SI
      bSI = √mu0SI     # magnetic unit in SI
      jSI = 1/√(mu0SI)   # current unit in SI
      c0 = 1.0             # speed of light (for Boris correction)
   elseif typeunit == "PLANETARY"
      xSI = 6378000         # Earth radius [default planet]
      tSI = 1.0             # s
      rhoSI = mpSI*1e6        # mp/cm^3
      uSI = 1e3             # km/s
      pSI = 1e-9            # nPa
      bSI = 1e-9            # nT
      jSI = 1e-6            # muA/m^2
      c0 = cSI/uSI         # speed of light in velocity units
   elseif typeunit == "OUTERHELIO"
      xSI = AuSI            # AU
      tSI = 1.0             # s
      rhoSI = mpSI*1e6        # mp/cm^3
      uSI = 1e3             # km/s
      pSI = 1e-1            # dyne/cm^2
      bSI = 1e-9            # nT
      jSI = 1e-6            # muA/m^2
      c0 = cSI/uSI         # speed of light in velocity units
   elseif typeunit == "SOLAR"
      xSI = RsSI            # radius of the Sun
      tSI = 1.0             # s
      rhoSI = 1e3             # g/cm^3
      uSI = 1e3             # km/s
      pSI = 1e-1            # dyne/cm^2
      bSI = 1e-4            # G
      jSI = 1e-6            # muA/m^2
      c0 = cSI/uSI         # speed of light in velocity units
   else
      throw(ArgumentError("invalid typeunit=$(typeunit)"))
   end

   # Overwrite values if given by eqpar
   for ieqpar in 1:neqpar
      var = param[ieqpar]
      if var == "xSI"
         xSI = eqpar[ieqpar]
      elseif var == "tSI"
         tSI = eqpar[ieqpar]
      elseif var == "uSI"
         uSI = eqpar[ieqpar]
      elseif var == "rhoSI"
         rhoSI = eqpar[ieqpar]
      elseif var == "mi"
         mi = eqpar[ieqpar]
      elseif var == "m1"
         m1 = eqpar[ieqpar]
      elseif var == "me"
         me = eqpar[ieqpar]
      elseif var == "qi"
         qi = eqpar[ieqpar]
      elseif var == "q1"
         q1 = eqpar[ieqpar]
      elseif var == "qe"
         qe = eqpar[ieqpar]
      elseif var == "g"
         gamma = eqpar[ieqpar]
      elseif var == "g1"
         gamma = eqpar[ieqpar]
      elseif var == "ge"
         ge = eqpar[ieqpar]
      elseif var == "c"
         c = eqpar[ieqpar]
      elseif var == "clight"
         clight = eqpar[ieqpar]
      elseif var == "r"
         r = eqpar[ieqpar]
      elseif var == "rbody"
         rbody = eqpar[ieqpar]
      end
   end

   # Overwrite distance unit if given as an argument
   if !isempty(distance)
      xSI = distance
   end

   # Overwrite ion & electron masses if given as an argument
   if !isempty(mp)
      Mi = mp
   end
   if !isempty(me)
      Me = me
   end

   # Calculate convenient conversion factors
   if typeunit == "NORMALIZED"
      ti0 = 1.0/Mi            # T      = p/rho*Mi           = ti0*p/rho
      cs0 = 1.0               # cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
      mu0A = 1.0               # vA     = sqrt(b/rho)        = sqrt(bb/mu0A/rho)
      mu0 = 1.0               # beta   = p/(bb/2)           = p/(bb/(2*mu0))
      uH0 = Mi                # uH     = j/rho*Mi           = uH0*j/rho
      op0 = 1.0/Mi            # omegap = sqrt(rho)/Mi       = op0*sqrt(rho)
      oc0 = 1.0/Mi            # omegac = b/Mi               = oc0*b
      rg0 = √Mi               # rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
      di0 = c0*Mi             # di = c0/sqrt(rho)*Mi        = di0/sqrt(rho)
      ld0 = Mi                # ld = sqrt(p)/(rho*c0)*Mi    = ld0*sqrt(p)/rho
   elseif typeunit == "PIC"
      ti0 = 1.0/Mi            # T      = p/rho*Mi           = ti0*p/rho
      cs0 = 1.0               # cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
      mu0A = 4*pi              # vA     = sqrt(b/(4*!pi*rho))= sqrt(bb/mu0A/rho)
      mu0 = 4*pi              # beta   = p/(bb/(8*!pi))     = p/(bb/(2*mu0))
      uH0 = Mi                # uH     = j/rho*Mi           = uH0*j/rho
      op0 = √(4π)/Mi          # omegap = sqrt(4*!pi*rho)/Mi = op0*sqrt(rho)
      oc0 = 1.0/Mi            # omegac = b/Mi               = oc0*b
      rg0 = √Mi               # rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
      di0 = 1.0/√(4π)         # di = 1/sqrt(4*!pi*rho)*Mi   = di0/sqrt(rho)
      ld0 = 1.0/√(4π)         # ld = sqrt(p/(4*!pi))/rho*Mi = ld0*sqrt(p)/rho
   else
      qom = eSI/(Mi*mpSI);
      moq = 1/qom
      ti0 = mpSI/kbSI*pSI/rhoSI*Mi       # T[K]=p/(nk) = ti0*p/rho
      cs0 = pSI/rhoSI/uSI^2              # cs          = sqrt(gs*p/rho)
      mu0A = uSI^2*mu0SI*rhoSI*bSI^(-2)   # vA          = sqrt(bb/(mu0A*rho))
      mu0 = mu0SI*pSI*bSI^(-2)           # beta        = p/(bb/(2*mu0))
      uH0 = moq*jSI/rhoSI/uSI            # uH=j/(ne)   = uH0*j/rho
      op0 = qom*√(rhoSI/e0SI)*tSI        # omegap      = op0*sqrt(rho)
      oc0 = qom*bSI*tSI                  # omegac      = oc0*b
      rg0 = moq*√(pSI/rhoSI)/bSI/xSI/√(Mi)    # rg     = rg0*sqrt(p/rho)/b
      di0 = cSI/(op0/tSI)/xSI                 # di=c/omegap = di0/sqrt(rho)
      ld0 = moq*√(pSI)/rhoSI/xSI              # ld          = ld0*sqrt(p)/rho
   end

   return true
end

function Base.show(io::IO, data::BATS)
   showhead(io, data)
   if data.list.bytes ≥ 1e9
      str = @sprintf "filesize: %.1f GB" data.list.bytes/1e9
   elseif data.list.bytes ≥ 1e6
      str = @sprintf "filesize: %.1f MB" data.list.bytes/1e6
   elseif data.list.bytes ≥ 1e3
      str = @sprintf "filesize: %.1f KB" data.list.bytes/1e3
   else
      str = @sprintf "filesize: %.1f Bytes" data.list.bytes
   end
   println(io, str)
   println(io, "snapshots: ", data.list.npictinfiles)
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
showhead(data::BATS) = showhead(data.list, data.head)
showhead(io::IO, data::BATS) = showhead(data.list, data.head, io)
