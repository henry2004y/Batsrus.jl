# Frequency analysis of 1D data.
#
# MAT package seems not able to read timeseries data. I can only work around
# this by identifying the starting and ending array index in magnetopause
# crossings.
#
# Hongyang Zhou, hyzhou@umich.edu 12/10/2019

using MAT, FFTW, Statistics, PyPlot

"""Moving box average."""
function sma(a::Array, n::Int)
   vals = zeros(size(a,1) - (n-1), size(a,2))

   for i in 1:size(a,1) - (n-1)
      for j in 1:size(a,2)
         vals[i,j] = mean(a[i:i+(n-1),j])
      end
   end

   vals
end

function B_FFT(filename, varname="Bobs", vartype="total", Fs=1/0.3,
   Brange=1500:6000)

   data = matread(filename)
   B = data[varname]
   if vartype == "total"
      Btmp = @. sqrt(B[Brange,1]^2 + B[Brange,2]^2 + B[Brange,3]^2)
   elseif vartype == "x"
      Btmp = B[Brange,1]
   elseif vartype == "y"
      Btmp = B[Brange,2]
   elseif vartype == "z"
      Btmp = B[Brange,3]
   end

   L = length(Btmp)

   Y_B = fft(Btmp)

   P2 = @. abs(Y_B/L)
   P1 = P2[1:floor(Int,L/2+1)]
   P1[2:end-1] .= 2 .* P1[2:end-1]

   f = Fs*(0:(L/2))/L

   return f, P1
end


filenames = ["/Users/hyzhou/Documents/research/Ganymede/data/B_G8.mat",
   "/Users/hyzhou/Documents/research/Ganymede/data/B_G8_PIC.mat",
   "/Users/hyzhou/Documents/research/Ganymede/data/B_G8_PIC_AMR2.mat",
   "/Users/hyzhou/Documents/research/Ganymede/data/B_G8_Hall_AMR2.mat"]

cType = "total" # ["x", "y", "z", "total"]

fObs, pObs   = B_FFT(filenames[1], "Bobs", cType)
fHall, pHall = B_FFT(filenames[1], "Bsim", cType, 1.0, 1573:2959)
fHall2, pHall2 = B_FFT(filenames[4], "Bsim", cType, 1.0, 1573:2959)
fPIC, pPIC   = B_FFT(filenames[2], "Bsim", cType, 1.0, 1573:2959)
fPIC2, pPIC2 = B_FFT(filenames[3], "Bsim", cType, 1.0, 1573:2959)

#=
figure()
plot(fObs[2:end], pObs[2:end].^2)
plot(fHall[2:end], pHall[2:end].^2)
plot(fPIC[2:end], pPIC[2:end].^2)
plt.yscale("log")
xlim(0.0, 0.3)
legend(["Obs","Hall","PIC"])
=#

nSmooth = 19
nSHalf = convert(Int,(nSmooth-1)/2)

figure()
plot(fObs[2+nSHalf:end-nSHalf], sma(pObs[2:end].^2, nSmooth),"k")
plot(fHall[2+nSHalf:end-nSHalf], sma(pHall[2:end].^2, nSmooth))
plot(fHall2[2+nSHalf:end-nSHalf], sma(pHall2[2:end].^2, nSmooth))
plot(fPIC[2+nSHalf:end-nSHalf], sma(pPIC[2:end].^2, nSmooth),"#d62728")
plot(fPIC2[2+nSHalf:end-nSHalf], sma(pPIC2[2:end].^2, nSmooth))
#plt.xscale("log")
plt.yscale("log")
xlim(0.0, 0.3)
grid(true)
legend(["Obs","Hall fine", "Hall coarse", "PIC fine","PIC coarse"])
