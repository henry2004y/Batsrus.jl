module VisAna
# Reimplementation of BATSRUS data reader in Julia
#
# Hongyang Zhou, hyzhou@umich.edu 07/24/2019

export readdata, plotdata, plotlogdata

using Glob
using PyPlot

struct Data
   x::Array{Float64}
   w::Array{Float64}
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
Filenames can be provided with wildcards. You are also allowed to pass multiple
filenames as a single string.
`fileheads, data, filelist = readdata(filename, "npict", 2, "verbose", false)`

# Examples
```jldoctest
filenames = "1d__raw*"
fileheads, data, filelist = BATSRUS.readdata(filenames)
```
"""
function readdata( filenamesIn::String; dir::String=".", npict::Int=1,
   verbose::Bool=true )

   ## Check the existence of files
   filenamesSplit = split(filenamesIn)
   filenames = Array{String,1}(undef,0)
   for filename in filenamesSplit
      filesfound = glob(filename)
      if isempty(filesfound)
         error("Error in readdata: no matching filename was found for
         $(filename)")
      end
      filenames = vcat(filenames, filesfound)
   end

   nfile = length(filenames)

   fileheads = []
   data = Vector{Data}()

   filelist, fileID, pictsize = get_file_types(nfile,filenames,dir)

   if verbose
      [println("filename=$(filelist[i].name)\n"*
      "npict=$(filelist[i].npictinfiles)") for i in length(filelist)]
   end

   for ifile=1:nfile
      if any(filelist[ifile].npictinfiles .- npict < 0)
         error("file $(ifile): npict out of range!")
      end
      seekstart(fileID[ifile])
   end

   ## Read data from files
   for ifile=1:nfile
      if occursin("log",filelist[ifile].type)
         filehead,data = readlogdata(joinpath(dir,filelist[ifile].name))
         return
      else
         # Skip npict-1 snapshots (because we only want npict snapshot)
         skip(fileID[ifile], pictsize[ifile]*(npict-1))

         filehead = getfilehead(fileID[ifile], filelist[ifile].type, 2)
         push!(fileheads, filehead)

         # Read data
         fileType = lowercase(filelist[ifile].type)
         if fileType == "ascii"
            x,w = getpictascii(fileID[ifile], fileheads[ifile])
         elseif fileType == "binary"
            x,w = getpictbinary(fileID[ifile], fileheads[ifile])
         elseif fileType == "real4"
            x,w = getpictreal(fileID[ifile], fileheads[ifile])
         elseif fileType == "log"

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

         push!(data,Data(x,w))

         println("Finished reading $(filelist[ifile].name)")
      end
      close(fileID[ifile])
   end

   return fileheads, data, filelist

end

"""
   read_log_data(filename)
Read information from log file.
"""
function read_log_data( filename )

# This part can be achieved using some packages.

#=
data = importdata(filename)

head = Dict(:ndim => 3, :headline => "", :it => Int32(-1.0),
   :time => Float32(-1.0), :gencoord => false, :neqpar => Int32(0), :nw => 1,
   :nx => [Int32(0)], :variables => Array{String,1}(undef,1))

if isstruct(data)
   head[:headline]  = data.textdata{1}
   head[:variables] = data.colheaders
   head[:ndim]      = 1
   head[:it]        = 0
   head[:time]      = 0.0
   head[:gencoord]  = false
   head[:nx]        = 1
   head[:nw]        = numel(data.colheaders)
   head[:variables] = split(data.textdata)

   data = []
end

return head, data
=#
end

"""
   get_file_types(nfile, filenames, dir)
Get the type of files.
...
# Output arguments
- `filelist::FileList`: fulfilled file structs.
- `fileID::Vector{IOStream}`: file IOStream for accessing data.
- `pictsize::Int`: size (in bytes) of one snapshot.
...
"""
function get_file_types(nfile::Int,filenames::Array{String,1},dir::String)

   fileID   = Array{IOStream,1}(undef,nfile)
   pictsize = Array{Int64,1}(undef,nfile)

   filelist = Array{FileList,1}(undef,nfile)

   for ifile=1:nfile
      f = joinpath(dir,filenames[ifile])
      fileID[ifile] = open(f,"r")

      bytes = filesize(filenames[ifile])
      type  = ""

      # Check the appendix of file names
      # I realized that this may not be a robust way; especially in linux
      # because you can have any appendix you want!!!
      # Gabor uses a trick: the first 4 bytes decides the file type()

      if occursin(r"^.*\.(log|dat)$",filenames[ifile])
         filelist[ifile][:type] = "log"
         filelist[ifile].npictinfiles = 1
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
               error("Error in get_file_types: strange unformatted file:
                  $(filelist[ifile].name)")
            end

            if lenhead==500
               type = uppercase(type)
            end
         end
      end

      # Obtain file size & number of snapshots
      seekstart(fileID[ifile])
      pictsize[ifile] = getfilehead(fileID[ifile], type)
      npictinfiles = floor(Int, bytes / pictsize[ifile])

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
function getfilehead(fileID::IOStream, type::String, iargout::Int=1)

   pictsize = 0

   head = Dict(:ndim => 3, :headline => "", :it => Int32(-1.0),
      :time => Float32(-1.0),
      :gencoord => false, :neqpar => Int32(0), :nw => 1, :nx => [Int32(0)],
      :eqpar => [Float32(0.0)], :variables => Array{String,1}(undef,1),
      :wnames => Array{String,1}(undef,1))

   # Create a struct for substituting common block file_head
   if isempty(head)
      println("Do something here?")
   end

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

   if iargout == 1
      return pictsize
   elseif iargout == 2
      # Set variables array
      head[:variables] = split(varname)     # returns a string array
      return head
   else
      error("unknown iargout=$(iargout)!")
   end

end

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
   plotdata(data, filehead, (func="rho", plotmode="contbar",
      plotrange=[-Inf Inf -Inf Inf],
      plotinterval=0.1))
Plot the variable from SWMF output.

`plotdata(data, filehead, "p", plotmode="contbar")`

`plotdata(data, filehead, "p", plotmode=grid)`

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
- `cut`: (optional) select 2D cut plane from 3D outputs ["x","y","z"].
- `CutPlaneIndex`: (optional)
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
   plotinterval::Float64=0.1, multifigure::Bool=true, verbose::Bool=true)

   x,w = data.x, data.w
   plotmode = split(plotmode)
   vars     = split(func)
   ndim     = filehead[:ndim]

   # I should check the size info of x to determine the type of plot!
   if ndim == 3 && isempty(cut)

   end

   if verbose
      println("============ PLOTTING PARAMETERS ===============")
      println("wnames = $(filehead[:wnames])")
      println("================================================")
      # Display min & max for each variable
      for var in vars
         VarIndex_ = findfirst(x->x==var,filehead[:wnames])
         if ndim == 1
            println("Min & Max value for $(var) :$(minimum(w[:,VarIndex_]))",
               ", $(maximum(w[:,VarIndex_]))")
         elseif ndim == 2
            println("Min & Max value for $(var) :$(minimum(w[:,:,VarIndex_]))",
               ", $(maximum(w[:,:,VarIndex_]))")
         end
      end
   end

   ## plot multiple variables with same plotmode
   if length(plotmode) < length(vars)
      [push!(plotmode, plotmode[i]) for i = 1:length(vars)-length(plotmode)]
   end

   ## Plot
   if ndim == 1
      for (ivar,var) in enumerate(vars)
         VarIndex_ = findfirst(x->x==var,filehead[:wnames])
         multifigure && figure()
         if plotmode[ivar] !== "scatter"
            plot(x,w[:,VarIndex_])
         else
            scatter(x,w[:,VarIndex_])
         end
         xlabel("x")
         ylabel("$(var)")
      end
   elseif ndim == 2
      for (ivar,var) in enumerate(vars)
         # I need to think of a better way to check. now this cannot identify the
         # vars for streamline | quiver plotting!!!
         if plotmode[ivar] ∈ ("mesh","meshbar","meshbarlog","cont","contf",
            "contbar","contlog","contbarlog")

            VarIndex_ = findfirst(x->x==var,filehead[:wnames])
            isempty(VarIndex_) && error("$(var) not found in header variables!")

            if multifigure figure() end

            if filehead[:gencoord] # Generalized coordinates
               # Reorganize & pick data in plot region only
               W = reshape(w[:,:,VarIndex_],:,1)
               X = reshape(x[:,:,1],:,1)
               Y = reshape(x[:,:,2],:,1)

               if all(abs.(plotrange) .!== Inf)
                  axis(plotrange)
               else
                  plotrange[1] = minimum(X)
                  plotrange[2] = maximum(X)
                  plotrange[3] = minimum(Y)
                  plotrange[4] = maximum(Y)
               end

               xyIndex = X .> plotrange[1] & X .< plotrange[2] & Y .> plotrange[3] & Y .< plotrange[4]
               X = X[xyIndex]
               Y = Y[xyIndex]
               W = W[xyIndex]
               # Default is linear interpolation
               F = scatteredInterpolant[X,Y,W]
               xq, yq = meshgrid(plotrange[1]:plotinterval:plotrange[2],plotrange[3]:plotinterval:plotrange[4])
               vq = F(xq,yq)

            else # Cartesian coordinates
               if all(abs.(plotrange) .== Inf)
                  xq = x[:,:,1]
                  yq = x[:,:,2]
                  vq = w[:,:,VarIndex_]
               else
                  # Note: the number of points in each dimension can be
                  # changed to the proper number you like. Here I use the
                  # number of points in each dimension from filehead.
                  xlin = range(plotrange[1],plotrange[2],length=filehead.nx[1])
                  ylin = range(plotrange[3],plotrange[4],length=filehead.nx[2])
                  xq, yq = meshgrid(xlin,ylin)
                  # 3D? 2D? My understanding is that as long as the last index
                  # is 1; it is fine.
                  X = x[:,:,:,1]
                  Y = x[:,:,:,2]
                  # From ndgrid to meshgrid format()
                  X  = permute(X,[2 1 3])
                  Y  = permute(Y,[2 1 3])
                  W  = permute(w[:,:,VarIndex_],[2 1 3])
                  vq = interp2(X,Y,W,xq,yq)
               end
            end

            if plotmode[ivar] == "contbar"
               contourf(xq,yq,vq,20,"Edgecolor','none"); c = colorbar()
            elseif plotmode[ivar] == "cont"
               contour(xq,yq,vq,20)
            elseif plotmode[ivar] == "contf"
               contourf(xq,yq,vq,20,"Edgecolor','none")
            elseif plotmode[ivar] == "contlog"
               contourf(xq,yq,log10(vq),20,"Edgecolor','none")
            elseif plotmode[ivar] == "contbarlog"
               contourf(xq,yq,log10(vq),20,"Edgecolor','none")
               c = colorbar()
               #c.Label.String = "log10"
            elseif plotmode[ivar] == "meshbar"
               mesh(xq,yq,vq); colorbar()
            elseif plotmode[ivar] == "mesh"
               mesh(xq,yq,vq)
            elseif plotmode[ivar] == "meshbarlog"
               mesh(xq,yq,log10(vq)); c = colorbar()
               #c.Label.String = "log10"
            end

            xlabel(filehead[:variables][1]); ylabel(filehead[:variables][2])
            title(filehead[:wnames][VarIndex_])
            dim = [0.125, 0.013, 0.2, 0.045]
            #str = sprintf("it=$(filehead[:it]), time=$(filehead[:time])")
            #annotation["textbox',dim,'String',str,'FitBoxToText','on","FontWeight','bold"]

         elseif plotmode[ivar] ∈ ("trimesh","trisurf","tricont","tristream")
            # triangular mesh()
            figure()
            # find the index for var in filehead.wnames
            VarIndex_ = findfirst(x->x==vars[ivar],filehead[:wnames])
            X = reshape(x[:,:,1],[],1)
            Y = reshape(x[:,:,2],[],1)
            W = reshape(w[:,:,VarIndex_],[],1)
            if !isempty(plotrange)
               xyIndex = X .> plotrange[1] & X .< plotrange[2] & Y .> plotrange[3] & Y .< plotrange[4];
               X = X[xyIndex]
               Y = Y[xyIndex]
               W = W[xyIndex]

               #             #tristream test
               #             u = reshape(w[:,:,2],1,[])
               #             v = reshape(w[:,:,4],1,[])
               #             u = u[xyIndex]
               #             v = v[xyIndex]
            end
            t = delaunayn([X, Y])
            if plotmode[ivar] == "trimesh"
               trimesh(t,X,Y,W)
            elseif plotmode[ivar] == "trisurf"
               trisurf(t,X,Y,W,EdgeColor="none")
            elseif plotmode[ivar] == "tristream"
               #             x0 = ones(7,1)*-3
               #             y0 = [ -3 -2 -1 0 1 2 3]'
               #             T = triangulation[t,x,y]
               #             #triplot(T)
               #             #FlowP=TriStream[T,u,v,x0,y0,1,2e3]
               #             FlowP=TriStream[T,u,v,x0,y0]
               #             PlotTriStream[FlowP,'r']
            else
               # I need to use patch to write my own tricont function!
               #tricont[t,x,y,w]
            end

         elseif plotmode[ivar] ∈ ("stream","streamover")
            if plotmode[ivar] == "streamover"
               # Overplotting with more variables; keyword "over"
            else
               if multifigure figure() end
            end

            # find the index for var in filehead.wnames
            VarStream  = split[func[ivar],';']
            VarIndexS1 = strcmpi(VarStream[1],filehead.wnames)
            VarIndexS2 = strcmpi(VarStream[2],filehead.wnames)

            if filehead.gencoord # Generalized coordinates
               F1 = scatteredInterpolant(x[:,1,1],x[:,1,2],w[:,1,VarIndexS1])
               v1 = F1(xq,yq)
               F2 = scatteredInterpolant(x[:,1,1],x[:,1,2],w[:,1,VarIndexS2])
               v2 = F2(xq,yq)
            else # Cartesian coordinates
               # 3D? 2D?
               X = x[:,:,:,1]
               Y = x[:,:,:,2]
               F1 = griddedInterpolant(X,Y,w[:,:,1,VarIndexS1])
               v1 = F1(xq,yq)
               F2 = griddedInterpolant(X,Y,w[:,:,1,VarIndexS2])
               v2 = F2(xq,yq)
            end

            # Modify the density of streamlines if needed
            s = streamslice(xq,yq,v1,v2,streamdensity,"linear")

            for is=1:length(s)
               s[is].Color = "w"; # Change streamline color to white
               s[is].LineWidth = 1.5
            end

            if !any(abs.(plotrange) .== Inf) axis(plotrange) end

         elseif plotmode[ivar] ∈ ("quiver","quiverover")
            if plotmode[ivar] == "quiverover"
               # Overplotting with more variables; keyword "over"
            else
               if multifigure figure() end
            end

            # find the index for var in filehead.wnames
            VarQuiver  = split(func[ivar],';')
            VarIndexS1 = findfirst(x->x==VarQuiver[1], filehead[:wnames])
            VarIndexS2 = findfirst(x->x==VarQuiver[2], filehead[:wnames])

            q = quiver(x[:,1,1],x[:,1,2],w[:,1,VarIndexS1],w[:,1,VarIndexS2])
            #q.Color = 'w'
            q.AutoScaleFactor = 0.5
            if !any(abs.(plotrange) .== Inf) axis(plotrange) end

         elseif plotmode[ivar] == "grid"    # Grid plot()
            figure()
            scatter(x[:,1,1],x[:,1,2],".",LineWidth=0.1)

            if !any(abs.(plotrange) .== Inf) axis(plotrange) end

            xlabel(filehead.variables[1]); ylabel(filehead.variables[2])
            title("Grid illustration")
            dim = [0.125, 0.013, 0.1, 0.046]
            #str = sprintf("it=$(filehead[:it])")
            #annotation["textbox',dim,'String',str,'FitBoxToText','on"]
         else
            error("unknown plot mode: $(plotmode[ivar])")
         end

      end

   else # 2D cut from 3D output; now only for Cartesian output
      for ivar = 1:length(vars)
         X = permute(x[:,:,:,1],[2 1 3])
         Y = permute(x[:,:,:,2],[2 1 3])
         Z = permute(x[:,:,:,3],[2 1 3])

         VarIndex_ = findfirst(x->x==vars[ivar],filehead[:wnames])
         if VarIndex_ == 0
            error("$(func[ivar]) not found in output variables!")
         end

         W  = permute(w[:,:,:,VarIndex_],[2 1 3])

         if multifigure figure() end

         if cut == "x"
            cut1 = squeeze(X[:,CutPlaneIndex,:])
            cut2 = squeeze(Z[:,CutPlaneIndex,:])
            W    = squeeze(W[:,CutPlaneIndex,:])
         elseif cut ==  "y"
            cut1 = squeeze(X[CutPlaneIndex,:,:])
            cut2 = squeeze(Z[CutPlaneIndex,:,:])
            W    = squeeze(W[CutPlaneIndex,:,:])
         elseif cut == "z"
            cut1 = squeeze(X[:,:,CutPlaneIndex])
            cut2 = squeeze(Y[:,:,CutPlaneIndex])
            W    = squeeze(W[:,:,CutPlaneIndex])
         end
      end

      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)

      contourf(cut1,cut2,W)
      colorbar;# axis equal

      if cut == "x"
         xlabel("y"); ylabel("z")
      elseif cut == "y"
         xlabel("x"); ylabel("z")
      elseif cut =="z"
         xlabel("x"); ylabel("y")
      end

   end

end

"""
   subsurface(x, y, data, limits)
Extract subset of surface dataset.
This is a simplified version of subvolume.
"""
function subsurface(x, y, data::Data, limits)

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

   if isnan(limits[1])  limits[1] = minimum(hx) end
   if isnan(limits[3])  limits[3] = minimum(hy) end
   if isnan(limits[2])  limits[2] = maximum(hx) end
   if isnan(limits[4])  limits[4] = maximum(hy) end

   xind = findfirst(limits[1]<=hx & hx<=limits[2])
   yind = findfirst(limits[3]<=hy & hy<=limits[4])

   newdata = subdata(data, xind, yind, sz)

   newx = x[xind, yind]
   newy = y[xind, yind]

   return  newx, newy, newdata
end

"""
   subdata(data, xind, yind, sz)
"""
function subdata(data, xind, yind, sz)

   newdata = data[xind, yind]
   newsz = size(newdata)

   if length(sz) > 2
      newdata = reshape(newdata, [newsz[1:3] sz[4:end]])
   end

   return newdata
end

"""
   animatedata(data, filehead)
Generate animations from data.
"""
function animatedata(data, filehead)

end


end
