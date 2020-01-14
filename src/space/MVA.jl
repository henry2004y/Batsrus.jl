# Minimum Variance Analysis (MVA)
#
# If the λ1 > λ2 > λ3 are 3 eigenvalues of the matrix M, then a good indicator
# of nice fitting LMN coordinate system should have
# λ2/λ3 > 5
#
# Example from Galileo satellite Ganymede flyby mag data.
# G8 inbound magnetopausee LMN for interval 15:51:40-15:52:20
# I cannot repeat the results of Figure 6 in [Kivelson 1998]
# -0.72 -0.55 0.43 from [Kivelson 1998]
# -0.942 -0.187 0.278 from elliptic empirical model
#
# Hongyang Zhou, hyzhou@umich.edu 01/10/2020

using Statistics, LinearAlgebra, CSV, Dates

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

	# Compute the eigen values and ratios (descending order)
	F = eigen(M, sortby = x -> -abs(x))
end

DoPlot = false
# Load magnetometer data Bx, By, Bz
filename = "../../test/Galileo_G8_flyby_MAG.dat"
df = CSV.File(filename; header=2, delim=" ", ignorerepeated=true)

# inbound crossing interval
index_ = 2656:2775
F = MVA(df.Bx[index_], df.By[index_], df.Bz[index_])

t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )

println("Eigenvalues:",F.values)
println("Eigenvectors:")
println(F.vectors[:,1])
println(F.vectors[:,2])
println(F.vectors[:,3])
println("Ratio of intermediate variance to minimum variance = ",
	round(F.values[2]/F.values[3],digits=3))


if DoPlot
	using PyPlot
	# Original Cartesian coordinates
	figure(figsize=(12,4))
	subplot(411)
	plot(t, df.Bx, label=L"B_x")
	legend()
	subplot(412)
	plot(t, df.By, label=L"B_y")
	legend()
	subplot(413)
	plot(t, df.Bz, label=L"B_z")
	legend()
	subplot(414)
	plot(sqrt.(df.Bx.^2 .+ df.By.^2 .+ df.Bz.^2), label=L"B")
	legend()
	tight_layout()

	# Coordinate transformation
	id_ = 2000:2900
	BL = sum(F.vectors[:,1]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
	BM = sum(F.vectors[:,2]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
	BN = sum(F.vectors[:,3]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
	figure(figsize=(12,4))
	subplot(411)
	plot(t[id_], BL, label=L"B_L")
	legend()
	subplot(412)
	plot(t[id_], BM, label=L"B_M")
	legend()
	subplot(413)
	plot(t[id_], BN, label=L"B_N")
	legend()
	subplot(414)
	plot(t[id_], sqrt.(BL.^2 .+ BM.^2 .+ BN.^2), label=L"B")
	legend()
	tight_layout()
end
