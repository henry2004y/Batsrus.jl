# Python calls Julia
import julia
j = julia.Julia()
j.include('VisAna.jl')
readdata = j.eval('VisAna.readdata')
filename='1d_bin.out'
filehead, data, filelist = readdata(filename,verbose=False);

# It turns out that I was too optimistic about using Python here. Even is the
# function calls are working, the data was saved as jlwrap object, which I
# cannot do any operation on it in Python.
