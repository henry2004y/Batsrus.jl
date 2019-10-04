module VisAna
# Reimplementation of BATSRUS data reader in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu

export readdata, readlogdata, plotdata, plotlogdata, animatedata, readtecdata
export Data, FileList, convertVTK

using Glob, PyPlot, Printf, PyCall, Dierckx, WriteVTK


struct Data{T}
   x::Array{T}
   w::Array{T}
end

struct FileList
   name::String
   type::String
   bytes::Int64
   npictinfiles::Int64
end

"""
   readdata(filenames,(, dir=".", npict=1, verbose=true))

Read data from BATSRUS output files. Stores the npict-th snapshot from an ascii
or binary data file into the x [coordinates] and w [data] arrays.
Filenames can be provided with wildcards.
`fileheads, data, filelist = readdata(filename, npict=2, verbose=false)`

# Examples
```jldoctest
filenames = "1d__raw*"
fileheads, data, filelist = readdata(filenames)
```
"""
function readdata( filenamesIn::String; dir::String=".", npict::Int=1,
   verbose::Bool=true )

   ## Check the existence of files
   filenames = Vector{String}(undef,0)
   filesfound = glob(filenamesIn, dir)
   if isempty(filesfound)
      error("readdata: no matching filename was found for $(filenamesIn)")
   end
   filenames = vcat(filenames, filesfound)

   nfile = length(filenames)

   fileheads = Vector{Dict}(undef,0)
   data = Vector{Data}()

   filelist, fileID, pictsize = getFileTypes(nfile,filenames)

   if verbose
      [println("filename=$(filelist[i].name)\n"*
      "npict=$(filelist[i].npictinfiles)") for i in eachindex(filelist)]
   end

   for ifile=1:nfile
      if any(filelist[ifile].npictinfiles .- npict < 0)
         error("file $(ifile): npict out of range!")
      end
      seekstart(fileID[ifile])
   end

   ## Read data from files
   for ifile=1:nfile
      # Skip npict-1 snapshots (because we only want npict snapshot)
      skip(fileID[ifile], pictsize[ifile]*(npict-1))

      filehead = getfilehead(fileID[ifile], filelist[ifile].type)
      push!(fileheads, filehead)

      # Read data
      fileType = lowercase(filelist[ifile].type)
      if fileType == "ascii"
         x,w = getpictascii(fileID[ifile], fileheads[ifile])
      elseif fileType == "binary"
         x,w = getpictbinary(fileID[ifile], fileheads[ifile])
      elseif fileType == "real4"
         x,w = getpictreal(fileID[ifile], fileheads[ifile])
      else
         error("get_pict: unknown filetype: $(filelist[ifile].type)")
      end

      setunits(fileheads[ifile],"")

      if verbose
         showhead(filelist[ifile], ifile, fileheads[ifile])
      end

      # Produce a wnames from the last file
      fileheads[ifile][:wnames] =
      fileheads[ifile][:variables][
      fileheads[ifile][:ndim]+1:fileheads[ifile][:ndim]+
      fileheads[ifile][:nw] ]

      push!(data, Data(x,w))

      println("Finished reading $(filelist[ifile].name)")

      close(fileID[ifile])
   end

   return fileheads, data, filelist
end

"""
   readlogdata(filename)
Read information from log file.
"""
function readlogdata( filename::String )
   # Is this really necessary?
   head = Dict(:ndim => 3, :headline => "", :it => Int32(-1.0),
   :time => Float32(-1.0), :gencoord => false, :neqpar => Int32(0), :nw => 1,
   :nx => [Int32(0)], :variables => Array{String,1}(undef,1))

   f = open(filename, "r")
   nLine = countlines(f)-2
   seekstart(f)
   head[:headline]  = readline(f)
   head[:variables] = split(readline(f))
   head[:ndim]      = 1
   head[:it]        = 0
   head[:time]      = 0.0
   head[:gencoord]  = false
   head[:nx]        = 1
   head[:nw]        = length(head[:variables])

   data = zeros(head[:nw],nLine)
   for i = 1:nLine
      line = split(readline(f))
      data[:,i] = parse.(Float64,line)
   end

   close(f)

   return head, data
end

"""
   readtecdata(filename, IsBinary, verbose)

Return header, data and connectivity from BATSRUS Tecplot outputs. Both binary
and ascii formats are supported.
# Examples
```jldoctest
filenames = "3d_ascii.dat"
fileheads, data, filelist = readtecdata(filenames)
```
"""
function readtecdata(filename::String, IsBinary::Bool=false,
   verbose::Bool=false)

   head = Dict(:variables => Array{String,1}(undef,1),
      :nNode => 0, :nCell => 0, :nDim  => 0)

   f = open(filename)

   # Read Tecplot header
   ln = readline(f)
   if startswith(ln, "TITLE")
      title = ln[8:end-1]
	  try # This should throw error if in binary format
		 parse(Int32,ln[17])
	  catch
		 IsBinary = true
	  end
   else
      @warn "No title provided."
   end
   ln = readline(f)
   if startswith(ln, "VARIABLES")
      VARS = split(ln[12:end],", ")
      for i in 1:length(VARS)
         VARS[i] = chop(VARS[i], head = 1, tail = 1)
      end
   else
      @warn "No variable names provided."
   end
   ln = readline(f)
   if startswith(ln, "ZONE")
      info = split(ln[6:end],", ")
	  nDim = parse(Int,info[1][4])
	  if IsBinary
		 nNode = read(IOBuffer(info[2][3:end]), Int32)
	     nCell = read(IOBuffer(info[3][3:end]), Int32)
	  else
      	 nNode = parse(Int,info[2][3:end])
	  	 nCell = parse(Int,info[3][3:end])
	 end
   else
      @warn "No zone info provided."
   end

   pt0 = position(f)
   while true
      x = readline(f)
      if !startswith(x, "AUXDATA")
         break
      else
         verbose && println(x)
      end
      pt0 = position(f)
   end
   seek(f, pt0)

   if IsBinary
      data = Array{Float32,2}(undef,length(VARS),nNode)
   else
   	  data = Array{Float32,2}(undef,length(VARS),nNode)
   end
   if nDim == 3
	  connectivity = Array{Int32,2}(undef,8,nCell)
   elseif nDim == 2
	  connectivity = Array{Int32,2}(undef,4,nCell)
   end

   if IsBinary
	  @inbounds for i = 1:nNode
		 read!(f, view(data,:,i))
	  end
	  @inbounds for i = 1:nCell
         read!(f, view(connectivity,:,i))
      end
   else
   	  @inbounds for i = 1:nNode
         x = readline(f)
         data[:,i] .= parse.(Float32, split(x))
      end
	  @inbounds for i = 1:nCell
		 x = readline(f)
		 connectivity[:,i] .= parse.(Int32, split(x))
	  end
   end

   close(f)

   head[:variables] = VARS
   head[:nNode] = nNode
   head[:nCell] = nCell
   head[:nDim]  = nDim

   return head, data, connectivity
end


"""
   getFileTypes(nfile, filenames, dir)
Get the type of files.
...
# Output arguments
- `filelist::FileList`: fulfilled file structs.
- `fileID::Vector{IOStream}`: file IOStream for accessing data.
- `pictsize::Int`: size (in bytes) of one snapshot.
...
"""
function getFileTypes(nfile::Int, filenames::Array{String,1})

   fileID   = Vector{IOStream}(undef,nfile)
   pictsize = Vector{Int64}(undef,nfile)

   filelist = Vector{FileList}(undef,nfile)

   for ifile=1:nfile
      f = filenames[ifile]
      fileID[ifile] = open(f,"r")

      bytes = filesize(filenames[ifile])
      type  = ""

      # Check the appendix of file names
      # Gabor uses a trick: the first 4 bytes decides the file type()

      if occursin(r"^.*\.(log)$",filenames[ifile])
         type = "log"
         npictinfiles = 1
      elseif occursin(r"^.*\.(dat)$",filenames[ifile])
         # tecplot ascii format
         type = "dat"
         npictinfiles = 1
      else
         # Obtain filetype based on the length info in the first 4 bytes
         lenhead = read(fileID[ifile], Int32)

         if lenhead!=79 && lenhead!=500
            type = "ascii"
         else
            # The length of the 2nd line decides between real4 & binary
            # since it contains the time; which is real*8 | real*4
            skip(fileID[ifile],lenhead+4)
            len = read(fileID[ifile],Int32)
            if len == 20
               type = "real4"
            elseif len == 24
               type = "binary"
            else
               error("Error in getFileTypes: strange unformatted file:
                  $(filelist[ifile].name)")
            end

            if lenhead==500
               type = uppercase(type)
            end
         end
         # Obtain file size & number of snapshots
         seekstart(fileID[ifile])
         pictsize[ifile] = getfilesize(fileID[ifile], type)
         npictinfiles = floor(Int, bytes / pictsize[ifile])
      end

      filelist[ifile] = FileList(filenames[ifile], type, bytes, npictinfiles)
   end

   return filelist, fileID, pictsize
end

"""
   getfilehead(fileID, type, iargout=1)
Obtain the header information from BATSRUS output files.
...
# Input arguments
- `fileID::IOStream`: file identifier.
- `type::String`: file type in ["ascii", "real4", "binary", "log"].
- `iargout::Int`: 1 for output pictsize, 2 for output filehead.
# Output arguments
- `pictsize::Int`: size of a single snapshot in bytes.
- `filehead::Dict`: file header info.
...
"""
function getfilehead(fileID::IOStream, type::String)

   pictsize = 0

   # Create a struct for substituting common block file_head
   head = Dict(:ndim => 3, :headline => "", :it => Int32(-1.0),
      :time => Float32(-1.0),
      :gencoord => false, :neqpar => Int32(0), :nw => 1, :nx => [Int32(0)],
      :eqpar => [Float32(0.0)], :variables => Array{String,1}(undef,1),
      :wnames => Array{String,1}(undef,1))

   ftype = string(lowercase(type))

   if ftype == type lenstr = 79 else lenstr = 500 end

   # Read header
   pointer0 = position(fileID)

   if ftype == "log"
      # This part is done in the main readdata function.
   elseif ftype == "ascii"
      head[:headline] = readline(fileID)
      line = readline(fileID)
      line = split(line)
      head[:it] = parse(Int,line[1])
      head[:time] = parse(Float64,line[2])
      head[:ndim] = parse(Int8,line[3])
      head[:neqpar] = parse(Int32,line[4])
      head[:nw] = parse(Int8,line[5])
      head[:gencoord] = head[:ndim] < 0
      head[:ndim] = abs(head[:ndim])
      line = parse(Int64, readline(fileID))
      head[:nx] = line
      if head[:neqpar] > 0
         head[:eqpar] = parse.(Float64, split(readline(fileID)))
      end
      varname = readline(fileID)

   elseif ftype ∈ ["real4","binary"]
      skip(fileID,4) # skip record start tag.
      head[:headline] = String(read(fileID, lenstr))
      skip(fileID,8) # skip record end/start tags.
      head[:it] = read(fileID,Int32)
      head[:time] = read(fileID,Float32)
      head[:ndim] = read(fileID,Int32)
      head[:gencoord] = (head[:ndim] .< 0)
      head[:ndim] = abs(head[:ndim])
      head[:neqpar] = read(fileID,Int32)
      head[:nw] = read(fileID,Int32)
      skip(fileID,8) # skip record end/start tags.
      head[:nx] = zeros(Int32,head[:ndim])
      read!(fileID,head[:nx])
      skip(fileID,8) # skip record end/start tags.
      if head[:neqpar] > 0
         head[:eqpar] = zeros(Float32,head[:neqpar])
         read!(fileID,head[:eqpar])
         skip(fileID,8) # skip record end/start tags.
      end
      varname = String(read(fileID, lenstr))
      skip(fileID,4) # skip record end tag.
   end

   # Header length
   pointer1 = position(fileID)
   headlen = pointer1 - pointer0

   # Calculate the snapshot size = header + data + recordmarks
   nxs = prod(head[:nx])

   if ftype == "log"
      pictsize = 1
   elseif ftype == "ascii"
      pictsize = headlen + (18*(head[:ndim]+head[:nw])+1)*nxs
   elseif ftype == "binary"
      pictsize = headlen + 8*(1+head[:nw]) + 8*(head[:ndim]+head[:nw])*nxs
   elseif ftype == "real4"
      pictsize = headlen + 8*(1+head[:nw]) + 4*(head[:ndim]+head[:nw])*nxs
   end

   # Set variables array
   head[:variables] = split(varname)     # returns a string array
   return head
end

""" Return the size of file. """
function getfilesize(fileID::IOStream, type::String)

   pictsize = 0

   ndim = convert(Int32,1)
   tmp  = convert(Int32,1)
   nw   = convert(Int32,1)
   nxs  = 0

   ftype = string(lowercase(type))

   if ftype == type lenstr = 79 else lenstr = 500 end

   # Read header
   pointer0 = position(fileID)

   if ftype == "log"
      # This part is done in the main readdata function.
   elseif ftype == "ascii"
      readline(fileID)
      readline(fileID)
      line = split(line)
      ndim = parse(Int32,line[3])
      tmp = parse(Int32,line[4])
      tmp > 0 && readline(fileID)
      nw = parse(Int32,line[5])
      readline(fileID)
      nx = parse(Int64, readline(fileID))
   elseif ftype ∈ ["real4","binary"]
      skip(fileID,4)
      read(fileID,lenstr)
      skip(fileID,8)
      read(fileID,Int32)
      read(fileID,Float32)
      ndim = abs(read(fileID,Int32))
      tmp = read(fileID,Int32)
      nw = read(fileID,Int32)
      skip(fileID,8)
      nx = zeros(Int32,ndim)
      read!(fileID,nx)
      skip(fileID,8)
      if tmp > 0
         tmp2 = zeros(Float32,tmp)
         read!(fileID,tmp2)
         skip(fileID,8) # skip record end/start tags.
      end
      read(fileID, lenstr)
      skip(fileID,4)
   end

   # Header length
   pointer1 = position(fileID)
   headlen = pointer1 - pointer0

   # Calculate the snapshot size = header + data + recordmarks
   nxs = prod(nx)

   if ftype == "log"
      pictsize = 1
   elseif ftype == "ascii"
      pictsize = headlen + (18*(ndim+nw)+1)*nxs
   elseif ftype == "binary"
      pictsize = headlen + 8*(1+nw) + 8*(ndim+nw)*nxs
   elseif ftype == "real4"
      pictsize = headlen + 8*(1+nw) + 4*(ndim+nw)*nxs
   end

   return pictsize
end

# There are plan to include this into Julia's base. See github for more info.
import Base: read!

"""
   read!(s,a)
Read slices of arrays using subarrays, in addition to the built-in methods.
"""
function read!(s::IO, a::SubArray{T}) where T

   for i in eachindex(a)
      a[i] = read(s, T)
   end
   return a
end

"""
   getpictascii(fileID, filehead)
Read ascii format data.
"""
function getpictascii(fileID::IOStream, filehead::Dict)

   ndim = filehead[:ndim]
   nw   = filehead[:nw]

   # Read coordinates & values row by row
   if ndim == 1 # 1D
      n1 = filehead[:nx][1]
      x  = Array{Float64,2}(undef,n1,ndim)
      w  = Array{Float64,2}(undef,n1,nw)
      for ix=1:n1
         temp = parse.(Float64, split(readline(fileID)))
         x[ix,:] .= temp[1]
         w[ix,:] .= temp[2:end]
      end
   elseif ndim == 2 # 2D
      n1 = filehead[:nx][1]
      n2 = filehead[:nx][2]
      x  = Array{Float64,3}(undef,n1,n2,ndim)
      w  = Array{Float64,3}(undef,n1,n2,nw)
      for ix1=1:n1, ix2=1:n2
         temp = parse.(Float64, split(readline(fileID)))
         x[ix1,ix2,:] .= temp[1]
         w[ix1,ix2,:] .= temp[2:end]
      end
   elseif ndim == 3 # 3D
      n1 = filehead[:nx][1]
      n2 = filehead[:nx][2]
      n3 = filehead[:nx][3]
      x  = Array{Float64,4}(undef,n1,n2,n3,ndim)
      w  = Array{Float64,4}(undef,n1,n2,n3,nw)
      for ix1=1:n1, ix2=1:n2, ix3=1:n3
         temp = parse.(Float64, split(readline(fileID)))
         x[ix1,ix2,ix3,:] .= temp[1]
         w[ix1,ix2,ix3,:] .= temp[2:end]
      end
   end

   return x, w
end

"""
   getpictbinary(fileID, filehead)
Read binary format data.
"""
function getpictbinary(fileID::IOStream, filehead::Dict)

   ndim = filehead[:ndim]
   nw   = filehead[:nw]

   # Read coordinates & values
   if filehead[:ndim] == 1 # 1D
      n1 = filehead[:nx][1]
      x  = Array{Float64,2}(undef,n1,ndim)
      w  = Array{Float64,2}(undef,n1,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw=1:nw
         read!(fileID, view(w,:,iw))
         skip(fileID,8) # skip record end/start tags.
      end
   elseif filehead[:ndim] == 2 # 2D
      n1 = filehead[:nx][1]
      n2 = filehead[:nx][2]
      x  = Array{Float64,3}(undef,n1,n2,ndim)
      w  = Array{Float64,3}(undef,n1,n2,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw=1:nw
         read!(fileID, view(w,:,:,iw))
         skip(fileID,8) # skip record end/start tags.
      end
   elseif filehead[:ndim] == 3 # 3D
      n1 = filehead[:nx][1]
      n2 = filehead[:nx][2]
      n3 = filehead[:nx][3]
      x  = Array{Float64,4}(undef,n1,n2,n3,ndim)
      w  = Array{Float64,4}(undef,n1,n2,n3,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw=1:nw
         read!(fileID, view(w,:,:,:,iw))
         skip(fileID,8) # skip record end/start tags.
      end
   end

   return x,w
end

"""
   getpictreal(fileID, filehead)
Read real4 format data.
"""
function getpictreal(fileID::IOStream, filehead::Dict)

   ndim = filehead[:ndim]
   nw   = filehead[:nw]

   # Read coordinates & values
   if filehead[:ndim] == 1 # 1D
      n1 = filehead[:nx][1]
      x  = Array{Float32,2}(undef,n1,ndim)
      w  = Array{Float32,2}(undef,n1,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw=1:nw
         read!(fileID, view(w,:,iw))
         skip(fileID,8) # skip record end/start tags.
      end
   elseif filehead[:ndim] == 2 # 2D
      n1 = filehead[:nx][1]
      n2 = filehead[:nx][2]
      x  = Array{Float32,3}(undef,n1,n2,ndim)
      w  = Array{Float32,3}(undef,n1,n2,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw=1:nw
         read!(fileID, view(w,:,:,iw))
         skip(fileID,8) # skip record end/start tags.
      end
   elseif filehead[:ndim] == 3 # 3D
      n1 = filehead[:nx][1]
      n2 = filehead[:nx][2]
      n3 = filehead[:nx][3]
      x  = Array{Float32,4}(undef,n1,n2,n3,ndim)
      w  = Array{Float32,4}(undef,n1,n2,n3,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw=1:nw
         read!(fileID, view(w,:,:,:,iw))
         skip(fileID,8) # skip record end/start tags.
      end
   end

   return x,w
end

"""
   setunits(filehead, type, (distunit, Mion, Melectron))
Set the units for the output files。
If type is given as "SI", "CGS", "NORMALIZED", "PIC", "PLANETARY", "SOLAR", set
typeunit = type otherwise try to guess from the fileheader.
Based on typeunit set units for distance [xSI], time [tSI], density [rhoSI],
pressure [pSI], magnetic field [bSI] and current density [jSI] in SI units.
Distance unit [rplanet | rstar], ion & electron mass in amu can be set with
optional distunit, Mion and Melectron.

Also calculate convenient constants ti0, cs0 ... for typical formulas.
This function needs to be improved!
"""
function setunits( filehead::Dict,type::String; distunit::Float64=1.0,
   Mion::Float64=1.0, Melectron::Float64=1.0)

   # This is currently not used, so return here
   return

   ndim      = filehead[:ndim]
   headline  = filehead[:headline]
   neqpar    = filehead[:neqpar]
   nw        = filehead[:nw]
   eqpar     = filehead[:eqpar]
   variables = filehead[:variables]

   mu0SI = 4*pi*1e-7      # H/m
   cSI   = 2.9978e8       # speed of light, [m/s]
   mpSI  = 1.6726e-27     # kg
   eSI   = 1.602e-19      # elementary charge, [C]
   AuSI  = 149597870700   # m
   RsSI  = 6.957e8        # m
   Mi    = 1.0            # Ion mass, [amu]
   Me    = 1.0/1836.15    # Electron mass, [amu]
   gamma = 5/3            # Adiabatic index for first fluid
   gammae= 5/3            # Adiabatic index for electrons
   kbSI  = 1.38064852e-23 # Boltzmann constant, [m2 kg s-2 K-1]
   e0SI  = 8.8542e-12     # [F/m]

   # This part is used to guess the units.
   # To be honest; I don`t understand the logic here. For example
   # nPa & m/s may appear in the same headline?
   if type !== ""
      typeunit = uppercase(type)
   elseif occursin("PIC",filehead[:headline])
      typeunit = "PIC"
   elseif occursin(" AU ", filehead[:headline])
      typeunit = "OUTERHELIO"
   elseif occursin(r"(kg/m3)|(m/s)", filehead[:headline])
      typeunit = "SI"
   elseif occursin(r"(nPa)|( nT )", filehead[:headline])
      typeunit = "PLANETARY"
   elseif occursin(r"(dyne)|( G)", filehead[:headline])
      typeunit = "SOLAR"
   else
      typeunit = "NORMALIZED"
   end

   if typeunit == "SI"
      xSI   = 1.0             # m
      tSI   = 1.0             # s
      rhoSI = 1.0             # kg/m^3
      uSI   = 1.0             # m/s
      pSI   = 1.0             # Pa
      bSI   = 1.0             # T
      jSI   = 1.0             # A/m^2
   elseif typeunit == "CGS"
      xSI   = 0.01            # cm
      tSI   = 1.0             # s
      rhoSI = 1000.0          # g/cm^3
      uSI   = 0.01            # cm/s
      pSI   = 0.1             # dyne/cm^2
      bSI   = 1.0e-4          # G
      jSI   = 10*cSI          # Fr/s/cm^2
   elseif typeunit == "PIC"
      # Normalized PIC units
      xSI   = 1.0             # cm
      tSI   = 1.0             # s
      rhoSI = 1.0             # g/cm^3
      uSI   = 1.0             # cm/s
      pSI   = 1.0             # dyne/cm^2
      bSI   = 1.0             # G
      jSI   = 1.0             # Fr/s/cm^2
      c0    = 1.0             # speed of light always 1 for iPIC3D
   elseif typeunit == "NORMALIZED"
      xSI   = 1.0             # distance unit in SI
      tSI   = 1.0             # time unit in SI
      rhoSI = 1.0             # density unit in SI
      uSI   = 1.0             # velocity unit in SI
      pSI   = 1.0             # pressure unit in SI
      bSI   = sqrt(mu0SI)     # magnetic unit in SI
      jSI   = 1/sqrt(mu0SI)   # current unit in SI
      c0    = 1.0             # speed of light (for Boris correction)
   elseif typeunit == "PLANETARY"
      xSI   = 6378000         # Earth radius [default planet]
      tSI   = 1.0             # s
      rhoSI = mpSI*1e6        # mp/cm^3
      uSI   = 1e3             # km/s
      pSI   = 1e-9            # nPa
      bSI   = 1e-9            # nT
      jSI   = 1e-6            # muA/m^2
      c0    = cSI/uSI         # speed of light in velocity units
   elseif typeunit == "OUTERHELIO"
      xSI   = AuSI            # AU
      tSI   = 1.0             # s
      rhoSI = mpSI*1e6        # mp/cm^3
      uSI   = 1e3             # km/s
      pSI   = 1e-1            # dyne/cm^2
      bSI   = 1e-9            # nT
      jSI   = 1e-6            # muA/m^2
      c0    = cSI/uSI         # speed of light in velocity units
   elseif typeunit == "SOLAR"
      xSI   = RsSI            # radius of the Sun
      tSI   = 1.0             # s
      rhoSI = 1e3             # g/cm^3
      uSI   = 1e3             # km/s
      pSI   = 1e-1            # dyne/cm^2
      bSI   = 1e-4            # G
      jSI   = 1e-6            # muA/m^2
      c0    = cSI/uSI         # speed of light in velocity units
   else
      error("invalid typeunit=$(typeunit)")
   end

   # Overwrite values if given by eqpar
   for ieqpar = 1:neqpar
      var = variables[ndim+nw+ieqpar]
      if var == "xSI"
         xSI   = eqpar[ieqpar]
      elseif var == "tSI"
         tSI   = eqpar[ieqpar]
      elseif var == "uSI"
         uSI   = eqpar[ieqpar]
      elseif var == "rhoSI"
         rhoSI = eqpar[ieqpar]
      elseif var == "mi"
         mi    = eqpar[ieqpar]
      elseif var == "m1"
         m1    = eqpar[ieqpar]
      elseif var == "me"
         me    = eqpar[ieqpar]
      elseif var == "qi"
         qi    = eqpar[ieqpar]
      elseif var == "q1"
         q1    = eqpar[ieqpar]
      elseif var == "qe"
         qe    = eqpar[ieqpar]
      elseif var == "g"
         gamma = eqpar[ieqpar]
      elseif var == "g1"
         gamma = eqpar[ieqpar]
      elseif var == "ge"
         ge    = eqpar[ieqpar]
      elseif var == "c"
         c     = eqpar[ieqpar]
      elseif var == "clight"
         clight= eqpar[ieqpar]
      elseif var == "r"
         r     = eqpar[ieqpar]
      elseif var == "rbody"
         rbody = eqpar[ieqpar]
      end
   end

   # Overwrite distance unit if given as an argument
   if !isempty(distunit) xSI = distunit end

   # Overwrite ion & electron masses if given as an argument
   if !isempty(Mion)      Mi = Mion end
   if !isempty(Melectron) Me = Melectron end

   # Calculate convenient conversion factors
   if typeunit == "NORMALIZED"
      ti0  = 1.0/Mi            # T      = p/rho*Mi           = ti0*p/rho
      cs0  = 1.0               # cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
      mu0A = 1.0               # vA     = sqrt(b/rho)        = sqrt(bb/mu0A/rho)
      mu0  = 1.0               # beta   = p/(bb/2)           = p/(bb/(2*mu0))
      uH0  = Mi                # uH     = j/rho*Mi           = uH0*j/rho
      op0  = 1.0/Mi            # omegap = sqrt(rho)/Mi       = op0*sqrt(rho)
      oc0  = 1.0/Mi            # omegac = b/Mi               = oc0*b
      rg0  = sqrt(Mi)          # rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
      di0  = c0*Mi             # di = c0/sqrt(rho)*Mi        = di0/sqrt(rho)
      ld0  = Mi                # ld = sqrt(p)/(rho*c0)*Mi    = ld0*sqrt(p)/rho
   elseif typeunit == "PIC"
      ti0  = 1.0/Mi            # T      = p/rho*Mi           = ti0*p/rho
      cs0  = 1.0               # cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
      mu0A = 4*pi              # vA     = sqrt(b/(4*!pi*rho))= sqrt(bb/mu0A/rho)
      mu0  = 4*pi              # beta   = p/(bb/(8*!pi))     = p/(bb/(2*mu0))
      uH0  = Mi                # uH     = j/rho*Mi           = uH0*j/rho
      op0  = sqrt(4*pi)/Mi     # omegap = sqrt(4*!pi*rho)/Mi = op0*sqrt(rho)
      oc0  = 1.0/Mi            # omegac = b/Mi               = oc0*b
      rg0  = sqrt(Mi)          # rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
      di0  = 1.0/sqrt(4*pi)    # di = 1/sqrt(4*!pi*rho)*Mi   = di0/sqrt(rho)
      ld0  = 1.0/sqrt(4*pi)    # ld = sqrt(p/(4*!pi))/rho*Mi = ld0*sqrt(p)/rho
   else
      qom  = eSI/(Mi*mpSI); moq = 1/qom
      ti0  = mpSI/kbSI*pSI/rhoSI*Mi       # T[K]=p/(nk) = ti0*p/rho
      cs0  = pSI/rhoSI/uSI^2              # cs          = sqrt(gs*p/rho)
      mu0A = uSI^2*mu0SI*rhoSI*bSI^(-2)   # vA          = sqrt(bb/(mu0A*rho))
      mu0  = mu0SI*pSI*bSI^(-2)           # beta        = p/(bb/(2*mu0))
      uH0  = moq*jSI/rhoSI/uSI            # uH=j/(ne)   = uH0*j/rho
      op0  = qom*sqrt(rhoSI/e0SI)*tSI     # omegap      = op0*sqrt(rho)
      oc0  = qom*bSI*tSI                  # omegac      = oc0*b
      rg0  = moq*sqrt(pSI/rhoSI)/bSI/xSI/sqrt(Mi) # rg     = rg0*sqrt(p/rho)/b
      di0  = cSI/(op0/tSI)/xSI                    # di=c/omegap = di0/sqrt(rho)
      ld0  = moq*sqrt(pSI)/rhoSI/xSI              # ld          = ld0*sqrt(p)/rho
   end

end

"""
   showhead(file, ifile, filehead)
Displaying file header information.
"""
function showhead(file::FileList, ifile::Int, filehead::Dict)

   println("----------------------")
   println("ifile     = $(ifile)")
   println("filename  = $(file.name)")
   println("filetype  = $(file.type)")
   println("headline  = $(filehead[:headline])")
   println("it        = $(filehead[:it])")
   println("time      = $(filehead[:time])")
   println("gencoord  = $(filehead[:gencoord])")
   println("ndim      = $(filehead[:ndim])")
   println("neqpar    = $(filehead[:neqpar])")
   println("nw        = $(filehead[:nw])")
   println("nx        = $(filehead[:nx])")
   println("----------------------")

   if filehead[:neqpar] > 0
      println("parameters = $(filehead[:eqpar])")
      println("coord names= $(filehead[:variables][1:filehead[:ndim]])")
      println("var   names= $(filehead[:variables][filehead[:ndim]+1:filehead[:ndim]+filehead[:nw]])")
      println("param names= $(filehead[:variables][filehead[:ndim]+filehead[:nw]+1:end])")
      println("=======================")
   end

end

"""
   plotlogdata(data, filehead, vars, (plotmode="line", plotrange=[-Inf,Inf]))
Plot information from log file.
...
# Input arguments
- `data::Data`: original variable data.
- `filehead::Dict`: header information.
- `vars::String`: variables for plotting.
- `plotmode::String`: (optional) type of plotting ["line","scatter"].
- `plotrange::Vector`: (optional) range of plotting.
...
"""
function plotlogdata(data::Data, filehead::Dict, func::String;
   plotmode::String="line", plotrange::Vector=[-Inf,Inf] )

   # This is intended for plotting log data, not data!

   vars     = split(func)
   plotmode = split(plotmode)

   for (ivar, var) in enumerate(vars)
      # find the index for var in filehead.variables
      VarIndex_ = findfirst(x->x==var, filehead[:variables])

      isnothing(VarIndex_) &&
      error("unknown plotting variable $(func[ivar])!")
      figure()
      if plotmode[ivar] == "line"
         plot(data[:,1],data[:,VarIndex_])
      elseif plotmode[ivar] == "scatter"
         scatter(data[:,1],data[:,VarIndex_])
      else
         error("unknown plot mode for plotlogdata!")
      end
      xlabel(filehead[:variables][1])
      ylabel(filehead[:variables][VarIndex_])
      title("log file data")
   end

end


"""
   plotdata(data, filehead, func, (...))
Plot the variable from SWMF output.

`plotdata(data, filehead, "p", plotmode="contbar")`

`plotdata(data, filehead, "p", plotmode="grid")`

`plotdata(data, filehead, func, plotmode="trimesh",plotrange=plotrange,
   plotinterval=0.2)`

...
# Input arguments
- `data::Data`: original variable data.
- `filehead::Dict`: header information.
- `vars::String`: variables for plotting.
- `plotmode::String`: (optional) type of plotting ["cont","contbar"]...
- `plotrange::Vector`: (optional) range of plotting.
- `plotinterval`: (optional) interval for interpolation.
- `density`: (optional) density for streamlines.
- `cut`: (optional) select 2D cut plane from 3D outputs ["x","y","z"].
- `cutPlaneIndex`: (optional)
- `streamdensity`: (optional) streamline density.
- `multifigure`: (optional) 1 for multifigure display, 0 for subplots.
- `verbose`: (optional) display additional information.
...
Right now this can only deal with 2D plots or 3D cuts. Full 3D plots may be
supported in the future.
I want to make this function more powerful to include plotting derived
variables, but it may not seem to be easy!
"""
function plotdata(data::Data, filehead::Dict, func::String; cut::String="",
   plotmode::String="contbar", plotrange::Vector{Float64}=[-Inf,Inf,-Inf,Inf],
   plotinterval::Float64=0.1, density::Float64=1.0, cutPlaneIndex::Int=1,
   multifigure::Bool=true, getrangeOnly::Bool=false, verbose::Bool=true)

   x,w = data.x, data.w
   plotmode = split(plotmode)
   vars     = split(func)
   ndim     = filehead[:ndim]
   nvar     = length(vars)

   if verbose || getrangeOnly
      println("============ PLOTTING PARAMETERS ===============")
      println("wnames = $(filehead[:wnames])")
      println("================================================")
      wmin = Vector{Float64}(undef,nvar)
      wmax = Vector{Float64}(undef,nvar)
      # Display min & max for each variable
      for (ivar,var) in enumerate(vars)
         if occursin(";",var) continue end # skip the vars for streamline
         VarIndex_ = findfirst(x->x==lowercase(var),
            lowercase.(filehead[:wnames]))
         if ndim == 1
            wmin[ivar] = minimum(w[:,VarIndex_])
            wmax[ivar] = maximum(w[:,VarIndex_])
         elseif ndim == 2
            wmin[ivar] = minimum(w[:,:,VarIndex_])
            wmax[ivar] = maximum(w[:,:,VarIndex_])
         end
         println("Min & Max value for $(var) :$(wmin[ivar])",", $(wmax[ivar])")
      end
      if getrangeOnly return wmin, wmax end
   end

   ## plot multiple variables with same plotmode
   if length(plotmode) < nvar
      [push!(plotmode, plotmode[i]) for i = 1:nvar-length(plotmode)]
   end

   ## Plot
   if ndim == 1
      for (ivar,var) in enumerate(vars)
         VarIndex_ = findfirst(x->x==var,filehead[:wnames])
         if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end
         if !occursin("scatter",plotmode[ivar])
            plot(x,w[:,VarIndex_])
         else
            scatter(x,w[:,VarIndex_])
         end
         if occursin("grid",plotmode[ivar])
            grid(true)
         end
         xlabel("x"); ylabel("$(var)")
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" filehead[:it] filehead[:time]
         at = matplotlib.offsetbox.AnchoredText(str,
                    loc="lower left", prop=Dict("size"=>8), frameon=true,
                    bbox_to_anchor=(0., 1.),
                    bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
      end
   elseif ndim == 2
      for (ivar,var) in enumerate(vars)
         occursin("over", plotmode[ivar]) && (multifigure = false)
         if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end
         if !occursin(";",var)
            VarIndex_ = findfirst(x->x==lowercase(var),
               lowercase.(filehead[:wnames]))
               isempty(VarIndex_) &&
                  error("$(var) not found in header variables!")
         end

         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")

            if filehead[:gencoord] # Generalized coordinates
               X = vec(x[:,:,1])
               Y = vec(x[:,:,2])
               W = vec(w[:,:,VarIndex_])

               if any(abs.(plotrange) .== Inf)
                  if plotrange[1] == -Inf plotrange[1] = minimum(X) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
               end

               # Create grid values first.
               xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
               yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)
               # Perform linear interpolation of the data (x,y) on grid(xi,yi)
               triang = matplotlib.tri.Triangulation(X,Y)
               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               np = pyimport("numpy")
               Xi, Yi = np.meshgrid(xi, yi)
               wi = interpolator(Xi, Yi)
            else # Cartesian coordinates
               if all(isinf.(plotrange))
                  xi = x[:,:,1]
                  yi = x[:,:,2]
                  wi = w[:,:,VarIndex_]
               else
                  if plotrange[1] == -Inf plotrange[1] = minimum(X) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

                  X = x[:,1,1]
                  Y = x[1,:,2]
                  W = w[:,:,VarIndex_]

                  xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
                  yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

                  spline = Spline2D(X, Y, W)
                  Xi = [i for i in xi, j in yi]
                  Yi = [j for i in xi, j in yi]
                  wi = spline(Xi[:], Yi[:])
                  wi = reshape(wi, size(Xi))'
               end
            end

            # I may need to use pattern match instead for a more robust method!
            if plotmode[ivar] == "contbar"
               c = contourf(xi,yi,wi,50)
            elseif plotmode[ivar] == "cont"
               c = contour(xi,yi,wi)
            elseif plotmode[ivar] == "contlog"
               c = contour(xi,yi,wi)
            elseif plotmode[ivar] == "contbarlog"
               c = contourf(xi,yi,wi)
            elseif plotmode[ivar] == "surfbar"
               c = plot_surface(xi,yi,wi)
            elseif plotmode[ivar] == "surfbarlog"
               c = plot_surface(xi,yi,wi)
            end

            occursin("bar", plotmode[ivar]) && colorbar()
            occursin("log", plotmode[ivar]) &&
               ( c.locator = matplotlib.ticker.LogLocator() )
            title(filehead[:wnames][VarIndex_])

         elseif plotmode[ivar] ∈ ("trimesh","trisurf","tricont","tristream")
            X = vec(x[:,:,1])
            Y = vec(x[:,:,2])
            W = vec(w[:,:,VarIndex_])

            # This needs to be modified!!!
            if !all(isinf.(plotrange))
               xyIndex = X .> plotrange[1] .& X .< plotrange[2] .&
                  Y .> plotrange[3] .& Y .< plotrange[4]
               X = X[xyIndex]
               Y = Y[xyIndex]
               W = W[xyIndex]
            end

            if plotmode[ivar] == "trimesh"
               triang = matplotlib.tri.Triangulation(X, Y)
               c = ax.triplot(triang)
            elseif plotmode[ivar] == "trisurf"
               c = ax.plot_trisurf(X, Y, W)
            elseif plotmode[ivar] == "tricont"
               c = ax.tricontourf(X, Y, W)
               fig.colorbar(c,ax=ax)
            elseif plotmode[ivar] == "tristream"
               error("not yet implemented!")
            end

            title(filehead[:wnames][VarIndex_])

         elseif plotmode[ivar] ∈ ("stream","streamover")
            VarStream  = split(var,";")
            VarIndex1_ = findfirst(x->x==lowercase(VarStream[1]),
               lowercase.(filehead[:wnames]))
            VarIndex2_ = findfirst(x->x==lowercase(VarStream[2]),
               lowercase.(filehead[:wnames]))

            if filehead[:gencoord] # Generalized coordinates
	            X, Y= vec(x[:,:,1]), vec(x[:,:,2])
	            if any(isinf.(plotrange))
                  if plotrange[1] == -Inf plotrange[1] = minimum(X) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end
               end

               # Create grid values first.
               xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
               yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

               # The PyCall here can be potentially replaced with Spline2D.
               # Perform linear interpolation of the data (x,y) on grid(xi,yi)
               triang = matplotlib.tri.Triangulation(X,Y)
               np = pyimport("numpy")
               Xi, Yi = np.meshgrid(xi, yi)
               W = w[:,1,VarIndex1_]

               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               v1 = interpolator(Xi, Yi)

               W = w[:,1,VarIndex2_]
               interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
               v2 = interpolator(Xi, Yi)

            else # Cartesian coordinates
               X = x[:,1,1]
               Y = x[1,:,2]
               if all(isinf.(plotrange))
                  Xi, Yi = X, Y
                  v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'
               else
                  if plotrange[1] == -Inf plotrange[1] = minimum(X) end
                  if plotrange[2] ==  Inf plotrange[2] = maximum(X) end
                  if plotrange[3] == -Inf plotrange[3] = minimum(Y) end
                  if plotrange[4] ==  Inf plotrange[4] = maximum(Y) end

                  w1, w2 = w[:,:,VarIndex1_], w[:,:,VarIndex2_]

                  xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
                  yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

                  Xi = [i for i in xi, j in yi]
                  Yi = [j for i in xi, j in yi]

                  spline = Spline2D(X, Y, w1)
                  v1 = spline(Xi[:], Yi[:])
                  v1 = reshape(v1, size(Xi))'

                  spline = Spline2D(X, Y, w2)
                  v2 = spline(Xi[:], Yi[:])
                  v2 = reshape(v2, size(Xi))'
               end
            end

            s = streamplot(Xi,Yi,v1,v2,color="w",linewidth=1.0,density=density)

         elseif occursin("quiver", plotmode[ivar])
            VarQuiver  = split(var,";")
            VarIndex1_ = findfirst(x->x==lowercase(VarQuiver[1]),
               lowercase.(filehead[:wnames]))
            VarIndex2_ = findfirst(x->x==lowercase(VarQuiver[2]),
               lowercase.(filehead[:wnames]))

            X, Y = x[:,1,1], x[1,:,2]
            v1, v2 = w[:,:,VarIndex1_]', w[:,:,VarIndex2_]'

            q = quiver(X,Y,v1,v2,color="w")

         elseif occursin("grid", plotmode[ivar])
            # This does not take subdomain plot into account!
            X, Y = x[:,:,1], x[:,:,2]
            scatter(X,Y,marker=".",alpha=0.6)
            title("Grid illustration")
         else
            error("unknown plot mode: $(plotmode[ivar])")
         end

         xlabel(filehead[:variables][1]); ylabel(filehead[:variables][2])
         dim = [0.125, 0.013, 0.2, 0.045]
         str = @sprintf "it=%d, time=%4.2f" filehead[:it] filehead[:time]
         at = matplotlib.offsetbox.AnchoredText(str,
                    loc="lower left", prop=Dict("size"=>8), frameon=true,
                    bbox_to_anchor=(0., 1.),
                    bbox_transform=ax.transAxes)
         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
         ax.add_artist(at)
         # recover status
         occursin("over", plotmode[ivar]) && (multifigure = true)
      end

   else # 2D cut from 3D output; now only for Cartesian output
      X = @view x[:,:,:,1]
      Y = @view x[:,:,:,2]
      Z = @view x[:,:,:,3]
      for (ivar,var) in enumerate(vars)
         if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")
            VarIndex_ = findfirst(x->x==lowercase(var),
               lowercase.(filehead[:wnames]))
            isempty(VarIndex_) && error("$(var) not found in header variables!")

            if ivar == 1 || multifigure fig, ax = subplots() else ax = gca() end

            W = w[:,:,:,VarIndex_]

            if cut ∈ ("x","")
               cut1 = @view X[cutPlaneIndex,:,:]
               cut2 = @view Y[cutPlaneIndex,:,:]
               W    = @view W[cutPlaneIndex,:,:]
            elseif cut ==  "y"
               cut1 = @view X[:,cutPlaneIndex,:]
               cut2 = @view Z[:,cutPlaneIndex,:]
               W    = @view W[:,cutPlaneIndex,:]
            elseif cut == "z"
               cut1 = @view X[:,:,cutPlaneIndex]
               cut2 = @view Y[:,:,cutPlaneIndex]
               W    = @view W[:,:,cutPlaneIndex]
            end
         elseif plotmode[ivar] ∈ ("stream","streamover")
            VarStream  = split(var,";")
            VarIndex1_ = findfirst(x->x==lowercase(VarStream[1]),
               lowercase.(filehead[:wnames]))
            VarIndex2_ = findfirst(x->x==lowercase(VarStream[2]),
               lowercase.(filehead[:wnames]))
            (isempty(VarIndex1_) || isempty(VarIndex2_)) &&
               error("$(VarStream) not found in header variables!")

            v1 = @view w[:,:,:,VarIndex1_]
            v2 = @view w[:,:,:,VarIndex2_]

            if cut ∈ ("x","")
               cut1 = @view Y[cutPlaneIndex,:,:]
               cut2 = @view Z[cutPlaneIndex,:,:]
               v1   = v1[cutPlaneIndex,:,:]'
               v2   = v2[cutPlaneIndex,:,:]'
            elseif cut ==  "y"
               cut1 = @view X[:,cutPlaneIndex,:]
               cut2 = @view Z[:,cutPlaneIndex,:]
               v1   = v1[:,cutPlaneIndex,:]'
               v2   = v2[:,cutPlaneIndex,:]'
            elseif cut == "z"
               cut1 = @view X[:,:,cutPlaneIndex]
               cut2 = @view Y[:,:,cutPlaneIndex]
               v1   = v1[:,:,cutPlaneIndex]'
               v2   = v2[:,:,cutPlaneIndex]'
            end
            cut1, cut2 = cut1', cut2'
         end

         if !all(isinf.(plotrange))
            cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
         end

	     if plotmode[ivar] ∈ ("surf","surfbar","surfbarlog","cont","contbar",
            "contlog","contbarlog")
            c = ax.contourf(cut1,cut2,W)
            fig.colorbar(c,ax=ax)
            #ax.axis("equal")

	     elseif plotmode[ivar] ∈ ("stream","streamover")
	        # Surprisingly, some box outputs do not have equal spaces???
	    	  #xi = range(cut1[1,1], stop=cut1[1,end], length=size(cut1)[2])
           #yi = range(cut2[1,1], stop=cut2[end,1], length=size(cut2)[1])

           xi = range(cut1[1,1], stop=cut1[1,end],
	       	  step=(cut1[1,end]-cut1[1,1])/(size(cut1,2)-1))
           yi = range(cut2[1,1], stop=cut2[end,1],
	           step=(cut2[end,1]-cut2[1,1])/(size(cut2,1)-1))

           Xi = [i for j in yi, x in xi]
           Yi = [j for j in yi, x in xi]

           s = streamplot(Xi,Yi,v1,v2,
              color="w",linewidth=1.0,density=density)
	     end

        if cut == "x"
           xlabel("y"); ylabel("z")
        elseif cut == "y"
           xlabel("x"); ylabel("z")
        elseif cut == "z"
           xlabel("x"); ylabel("y")
        end

     end

   end

end

"""
   subsurface(x, y, data, limits)
Extract subset of 2D surface dataset.
This is a simplified version of subvolume.
"""
function subsurface(x::Array{Float64,2}, y::Array{Float64,2},
   data::Array{Float64,2}, limits::Vector{Float64})

   if length(limits)!=4
      error("Reduction must be [xmin xmax ymin ymax]")
   end

   if limits[1] > limits[2]
      error("subvolume:InvalidReductionXRange")
   end
   if limits[3] > limits[4]
      error("subvolume:InvalidReductionYRange")
   end

   sz = size(data)

   hx = x[:,1]
   hy = y[1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])

   newdata = subdata(data, xind, yind, sz)

   newx = x[xind, yind]
   newy = y[xind, yind]

   return newx, newy, newdata
end

"""
   subdata(data, xind, yind, sz)
Return the sliced data based on indexes.
"""
function subdata(data::Array{Float64,2},
   xind::Vector{Int64}, yind::Vector{Int64}, sz::Tuple{Int64,Int64})

   newdata = data[xind, yind]
   newsz = size(newdata)

   if length(sz) > 2
      newdata = reshape(newdata, (newsz[1:3]..., sz[4:end]))
   end

   return newdata
end

"""
   animatedata(filelist, func, (plotmode="contbar",
      plotrange=[-Inf Inf -Inf Inf],
      plotinterval=0.1))
Generate animations from data. This is basically calling plotdata function for
multiple snapshots. The main issue here is to determine the colorbar/axis range
in advance to avoid any jump in the movie.
"""
function animatedata(filelist::FileList,func::String;
   imin::Int=1, imax::Int=1, cut::String="",
   plotmode::String="contbar", plotrange::Vector{Float64}=[-Inf,Inf,-Inf,Inf],
   plotinterval::Float64=0.1, verbose::Bool=true)


end


function animate(i,filelist)
   clf()
   fhead, d, flist = readdata(filelist.name,verbose=false,npict=i+1)
   plotdata(d[1],fhead[1],"p",plotmode="contbar")

   return gca()
end

"""
   convertVTK(head, data, connectivity, filename)

Convert 3D unstructured Tecplot data to VTK. Note that if using voxel type data
in VTK, the connectivity sequence is different from Tecplot.
"""
function convertVTK(head::Dict, data::Array{Float32,2},
   connectivity::Array{Int32,2}, filename::String="3DBATSRUS")

   nVar = length(head[:variables])

   points = @view data[1:head[:nDim],:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:nDim] == 3
      # PLT to VTK index_ = [1 2 4 3 5 6 8 7]
      for i = 1:2
         connectivity = swaprows(connectivity, 4*i-1, 4*i)
      end
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:nDim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid(filename, points, cells)

   for ivar = head[:nDim]+1:nVar
      if occursin("_x",head[:variables][ivar]) # vector
         var1 = @view data[ivar,:]
         var2 = @view data[ivar+1,:]
         var3 = @view data[ivar+2,:]
         namevar = replace(head[:variables][ivar], "_x"=>"")
         vtk_point_data(vtkfile, (var1, var2, var3), namevar)
      elseif occursin(r"(_y|_z)",head[:variables][ivar])
         continue
      else
         var = @view data[ivar,:]
         vtk_point_data(vtkfile, var, head[:variables][ivar])
      end
   end

   outfiles = vtk_save(vtkfile)
end

function swaprows(X::Array{Int32,2}, i::Int64, j::Int64)
   m, n = size(X)
   if (1 <= i <= n) && (1 <= j <= n)
      for k = 1:n
        @inbounds X[i,k],X[j,k] = X[j,k],X[i,k]
      end
      return X
   else
      throw(BoundsError())
   end
end

end
