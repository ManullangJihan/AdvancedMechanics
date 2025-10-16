# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")


# ---
# Coupled Spring mass damper
using DifferentialEquations
using Plots

function coupled_spring_mass_damper!(du, u, p, t)
    m1, m2, L1, L2, k1, k2, c1, c2 = p
    
    du[1] = u[3]
    du[2] = u[4]
    
    du[3] = -c1 * u[3] - k1 * (u[1] - L1) + k2 * (u[2] - u[1] - L2) / m1
    du[4] = (-c2 * u[4] - k2 * (u[2] - u[1] - L2)) / m2
end

# mass
m1 = 5
m2 = 2

# length
l1 = 3
l2 = 2

# spring constant
k1 = 0.3
k2 = 0.3

# friction constant
c1 = 0.1
c2 = 0.1

# Initial Condition
x₁0 = 1
x₂0 = x₁0+l2+1
v₁0 = 0
v₂0 = 0

u0 = [x₁0, x₂0, v₁0, v₂0]
params = (m1, m2, l1, l2, k1, k2, c1, c2)

t_end = 10.0
tspan = (0.0, t_end)
t = tspan[1]:0.1:tspan[2]

# Solve the ODE
prob = ODEProblem(coupled_spring_mass_damper!, u0, tspan, params)
sol = Array(solve(prob, Tsit5(), saveat=t))
x1s, x2s = sol[1, :], sol[2, :]

rectangle(x, y, w, h) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

W = 2.0
H = 1.0

@gif for (x1, x2) in zip(x1s, x2s)
	p = plot(rectangle(x1[1] - 1.0, 0.0, W, H), xlims=(-l1, 15.0), ylims=(0.0, 6.0),
		legend=false)
	plot!(p, [-l1, x1[1]-1.0], [0.5, 0.5], linewidth=3.0)

	plot!(p, [x1[1]+1.0, x2[1]+l2], [0.5, 0.5], linewidth=3.0)
	plot!(p, rectangle(x2[1]+l2, 0.0, W, H))
end

