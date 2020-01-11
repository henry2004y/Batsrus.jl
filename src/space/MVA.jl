# Minimum Variance Analysis (MVA)
#
# Example from Galileo satellite Ganymede flyby mag data.
# G8 inbound magnetopausee LMN for interval 15:51:40-15:52:20
#
# Hongyang Zhou, hyzhou@umich.edu 01/10/2020

using Statistics, LinearAlgebra, CSV, PyPlot

function MVA(Bx, By, Bz)

	B̄1 = mean(Bx)
	B̄2 = mean(By)
	B̄3 = mean(Bz)
	B̄11= mean(Bx.*Bx) - B̄1*B̄1
	B̄22= mean(By.*By) - B̄2*B̄2
	B̄33= mean(Bz.*Bz) - B̄3*B̄3
	B̄12= mean(Bx.*By) - B̄1*B̄2
	B̄23= mean(By.*Bz) - B̄2*B̄3
	B̄31= mean(Bz.*Bx) - B̄3*B̄1
	# Construct the matrix
	M = [B̄11 B̄12 B̄31; B̄12 B̄22 B̄23; B̄31 B̄23 B̄33]

	# Compute the eigen values and ratios (ascending order)
	F = eigen(M)
end

# Load magnetometer data Bx, By, Bz
filename = "../../test/Galileo_G8_flyby_MAG.dat"
df = CSV.File(filename; header=2, delim=" ", ignorerepeated=true)

# inbound crossing interval
index_ = 2656:2775
#index_ = 2686:2750
F = MVA(df.Bx[index_], df.By[index_], df.Bz[index_])

figure()
subplot(411)
plot(df.Bx)
subplot(412)
plot(df.By)
subplot(413)
plot(df.Bz)
#subplot(414)
#plot(sqrt.(BL.^2 .+ BM.^2 .+ BN.^2))

# Coordinate transformation

F.values

id_ = 2000:2900
BL = sum(F.vectors[:,3]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
BM = sum(F.vectors[:,2]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
BN = sum(F.vectors[:,1]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)

figure()
subplot(411)
plot(BL)
subplot(412)
plot(BM)
subplot(413)
plot(BN)
subplot(414)
plot(sqrt.(BL.^2 .+ BM.^2 .+ BN.^2))
# How do I check the correctness?
# There seems to be errors!
# -0.72 -0.55 0.43 from [Kivelson 1998]
# -0.942 -0.187 0.278 from elliptic empirical model
