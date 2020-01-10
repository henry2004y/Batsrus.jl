using FFTW, PyPlot, DelimitedFiles
using DSP, Statistics

include("/Users/hyzhou/Documents/Computer/Julia/FFT/savitskyGolay.jl")

const Rg = 2634000.0 #[m]

# Background
# G8 B= -10 -6 -86 [nT], U=140 0 0 [km/s]
Uxbk = 140; Bybk = -6; Bzbk = -86
TiltedAngle = atan(6/86)

# Get upstream theoretical maximum potential drop
filename = "../data/boundary.txt"
bc = readdlm(filename, ',', header=true)
data = bc[1]

#Δy =  maximum(x) - minimum(x)
Δy = 6.0 # [Rg]

# Take clock angle in yz plane into consideration
Potential_bk = Δy * cos(TiltedAngle) * Rg *
   Uxbk*1e3 * abs(Bzbk*cos(TiltedAngle)+Bybk*sin(TiltedAngle))*1e-9 *
   1e-3 # [kV]

Fs = 1           # Sampling frequency
T = 1/Fs         # Sampling period
L = 1197         # Length of signal
t = (0:L-1)*T    # Time vector

data_hall = readdlm("../data/CPCP_hall.txt")

# find the missing data and pad to them!!!

fig = figure(figsize=(8.5,4.2))
Y = welch_pgram(data_hall[:,2], 180, 100; onesided=true, fs=Fs)
plot(1 ./ Y.freq[2:end], Y.power[2:end])
plt.yscale("log")
plt.xscale("log")

data_pic  = readdlm("../data/CPCP_pic.txt")

Y = welch_pgram(data_pic[:,2], 200, 100; onesided=true, fs=Fs)
plot(1 ./ Y.freq[2:end], Y.power[2:end])
xlabel("Periods [s]", fontsize=12)
ylabel("Energy", fontsize=12)
grid("on")
legend(("Hall","PIC"))


println("mean Hall CPCP: ", mean(data_hall[:,2]))
println("mean PIC CPCP: ", mean(data_pic[:,2]))
