# Test of VisAna
#
# Hongyang Zhou, hyzhou@umich.edu 07/30/2019

# For better usage in REPL, I need to provide a script or wrapper to the full
# functions. For example, I want something like
# readdata
# in IDL. To achieve that, probably I need a script to do so.

using Test

include("VisAna.jl")
using .VisAna

filename = "1d_bin.out";

@testset "some tests" begin
   @test filehead, data, filelist = readdata(filename);
end

# 1D
filename = "1d_bin.out";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"p",plotmode="line")
plotdata(data[1],filehead[1],"p",plotmode="linegrid")
animatedata(data[1],filehead[1],"p")

# 2D Cartesian (structured)
filename = "z=0_raw_1_t10.07620_n00000102.out";
filehead, data, filelist = readdata(filename,verbose=false);

plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar streamover")
plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar quiverover")
plotdata(data[1],filehead[1],"p bx;by",plotmode="contbar streamover",
   density=2.0)
plotdata(data[1],filehead[1],"p",plotmode="grid")
plotdata(data[1],filehead[1],"p",plotmode="contbar",
   plotrange=[-50., 50., -1., 1.])
plotdata(data[1],filehead[1],"p",plotmode="contbar")
plotdata(data[1],filehead[1],"p",plotmode="contbarlog")
plotdata(data[1],filehead[1],"p",plotmode="surfbar")

# 2D unstructured
filename = "y=0_unstructured.outs";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar")
plotdata(data[1],filehead[1],"rho",plotmode="trimesh")
plotdata(data[1],filehead[1],"rho",plotmode="tricont")

animatedata(data[1],filehead[1],filelist[1],"rho",plotmode="contbar")

# 2D structured spherical ???
filename = "y_structured.out"; #???
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"rho",plotmode="contbar")


# 3D box
filename = "box.outs";
filehead, data, filelist = readdata(filename,verbose=false);
plotdata(data[1],filehead[1],"bx",plotmode="contbar",cut="y")
plotdata(data[1],filehead[1],"bx",plotmode="contbar",cut="y",
   plotrange=[-1.4,-1.1,0.70,0.78])

# 3D structured spherical
filename = "3d_structured.out";
filehead, data, filelist = readdata(filename,verbose=false);

# log file
logfilename = "shocktube.log";
filehead, data = readlogdata(logfilename)
