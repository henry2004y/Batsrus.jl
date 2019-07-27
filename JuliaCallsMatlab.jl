# Calling Matlab functions in Julia.
#
# This requires that Matlab can be called from the command line.
#
# Hongyang Zhou, hyzhou@umich.edu 07/17/2019

using MATLAB

filename="/Users/hyzhou/SWMF/test/BATSRUS/run_test/RESULT/GM/"*
   "z=0_raw_1_t00.00000_n00000000.out";

filehead, data = mxcall(:read_data,2,filename);

data = data["file1"];

mat"plot_data($data,$filehead,'p','plotmode','contbar')"
