# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using DifferentialEquations
using Plots

function triple_pendulum!(du, u, p, t)
    θ₁, ω₁, θ₂, ω₂, θ₃, ω₃ = u
    l₁, l₂, l₃, m₁, m₂, m₃, g = p
    
    c1 = cos(θ₁ - θ₂)
    c2 = cos(θ₂ - θ₃)
    s1 = sin(θ₁ - θ₂)
    s2 = sin(θ₂ - θ₃)
    
    denom1 = m₁ * l₁ - m₂ * l₁ * c1^2
    denom2 = m₂ * l₂ - m₃ * l₂ * c2^2
    
    du[1] = ω₁
    du[2] = (m₂ * l₁ * c1 * s1 * ω₁^2 + m₂ * g * s1 * c1 + m₂ * l₂ * c2 * s2 * ω₂^2 + m₂ * g * s2 * c2 + m₃ * l₃ * s2 * ω₃^2 + m₃ * g * s2) / denom1
    du[3] = ω₂
    du[4] = (-m₂ * l₂ * c2 * s2 * ω₂^2 - m₂ * g * s2 * c2 - m₃ * l₃ * s2 * ω₃^2 - m₃ * g * s2) / denom2
    du[5] = ω₃
    du[6] = -m₃ * l₃ * s2 * c2 * ω₃^2 - m₃ * g * s2 * c2 / (m₂ * l₂)
end

# Set parameters
l₁ = 1.0
l₂ = 1.0
l₃ = 1.0
m₁ = 1.0
m₂ = 1.0
m₃ = 1.0
g = 9.81

# Initial conditions
u₀ = [π/2, 0, π/2, 0, π/2, 0]  # Starting from a configuration where all pendulums are at rest at 90 degrees

# Time span
tspan = (0.0, 10.0)

# Set up ODE problem
prob = ODEProblem(triple_pendulum!, u₀, tspan, [l₁, l₂, l₃, m₁, m₂, m₃, g])

# Solve ODE
sol = solve(prob, Vern7())

# Plot the results
plot(sol, vars=(1, 3, 5), xlabel="θ₁", ylabel="θ₂", zlabel="θ₃", legend=false)
