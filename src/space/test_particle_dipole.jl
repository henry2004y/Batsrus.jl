# Particle pusher in Earth's dipole field.
# I am having issue with nondimensionalization!
#
# In DifferentialEquations, there are so called Ensemble simulations for many
# particles in parallel. Make use of that when going to many test particles.

using PyPlot
#using Plots
using DifferentialEquations
using LinearAlgebra # vector norm function, cross product
using Statistics # mean function

# Physical constants in SI units
const c = 3e8 #[m/s]
const q = 1.60217662e-19 #[C]
const mᵢ = 1.6726219e-27 #[kg]
const mₑ = 9.10938356e-31 #[kg]
const Rₑ = 6.38e6 #[m]
const μ = 4*π*1e-7  # Vacuum permeability, [V*s/(A*m)]

"Momentum equation for updating particle motion from EM force."
function derivatives!(dy::Vector{Float64},y::Vector{Float64},p::Vector{Any},
   t::Float64)
   q, m, E, BMoment = p
   dy[1:3] = y[4:6]
   dy[4:6] = q/m*(E + y[4:6] × (bfield(y[1:3],BMoment)))
end

"Calculates the magnetic field from a dipole M at position r."
function bfield(rIn::Vector{Float64},M::Vector{Float64})
   x,y,z = rIn
   r = sqrt(x^2 + y^2 + z^2)
   Coef = μ/(4*π*r^5)

   B = [3*x^2-r^2 3*x*y 3*x*z;
        3*y*x 3*y^2-r^2 3*y*z;
        3*z*x 3*z*y 3*z^2-r^2]*
        M*Coef

   return B
end

"Convert from spherical coordinates vector to Cartesian vector."
function sph2cart(r,ϕ,θ)
   r*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end

"""
   fieldline(L,ϕ,nP)

Creates points on one field line of the magnetic field from a dipole.
In a centered dipole magnetic field model, the path along a given L shell can be
described as
r = L*cos²λ,
where r is the radial distance (in planetary radii) to a point on the line,
λ is its co-latitude, and L is the L-shell of interest.
"""
function fieldline(ϕ::Float64,L::Float64=2.5,nP::Int=100)

   xyz = [ sph2cart(L*sin(θ)^2,ϕ,θ) for θ in range(-π,stop=π,length=nP) ]
   x = Vector{Float64}(undef,length(xyz))
   y = Vector{Float64}(undef,length(xyz))
   z = Vector{Float64}(undef,length(xyz))

   for (i, pos) in enumerate(xyz)
      x[i],y[i],z[i] = [pos...]
   end

   (x,y,z)
end

E = [0.0, 0.0, 0.0] # Electric field, [V/m]
#BMoment = [0.0, 0.0, 0.305e-4*Rₑ^2]
BMoment = [0.0, 0.0, 7.94e22] # [[V*s/(A*m)]
Ek = 5e7 # [eV]

# Velocity, [m/s]
v₀ = sph2cart(c*sqrt(1-1/(1+Ek*q/(mᵢ*c^2))^2), 0.0, π/4)
# Position, [m]
r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)

p = [q, mᵢ, E, BMoment]
tspan = (0.0, 10.0)
prob = ODEProblem(derivatives!, [r₀..., v₀...], tspan, p)
sol = solve(prob)

x = getindex.(sol.u,1) / Rₑ
y = getindex.(sol.u,2) / Rₑ
z = getindex.(sol.u,3) / Rₑ
u = getindex.(sol.u,4)
v = getindex.(sol.u,5)
w = getindex.(sol.u,6)

#=
# Visualization with Plots
# Plots is still not ready for productive usage...
pyplot()
display(plot3d(x,y,z,xlabel = "x", ylabel="y", zlabel="z",
   title="50 MeV ion trajectory in Earth dipole field"))

for ϕ in range(0, stop=2*π, length=10)
   display(plot!(fieldline(ϕ)...,color=:red, alpha=0.3, legend=false))
end
=#

# Visualization
fig = plt.figure()
using3D()
ax = plt.axes(projection="3d")

plot(x,y,z)

#mean_B = mean(norm.([bfield(sol'[i,1:3],BMoment) for i in 1:size(sol,2)]))

for ϕ in range(0, stop=2*π, length=10)
   plot(fieldline(ϕ)...,color=:red, alpha=0.3)
end
xlabel("x")
ylabel("y")
zlabel("z")
title("50 MeV ion trajectory in Earth dipole field")
