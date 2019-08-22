# Python calls Julia
import julia
j = julia.Julia()
j.include('VisAna.jl')
readdata = j.eval('VisAna.readdata')
filename='1d_bin.out'
filehead, data, filelist = readdata(filename,verbose=False);
