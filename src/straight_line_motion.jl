# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using Plots
using LinearAlgebra
using DifferentialEquations

# ---
# Straight Line Motion with no Force
function eqn1(du, u, p, t)
    du[1] = u[2]
    du[2] = 0
end

u01 = [1.0, 10.0]
tspan = (0.0, 10.0)
prob = ODEProblem(eqn1, u01, tspan)
sol = solve(prob)
plot(sol, lab=["Position" "Velocity"])


# ---
# Straight Line Motion with Constant Force
function eqn2(du, u, p, t)
    m = p[1]
    du[1] = u[2]
    du[2] = 10 / m # Constant Acceleration
end

mass2 = 5.0
p2 = [mass2]
u02 = [1.0, 0.0]
tspan2 = (0.0, 10.0)
prob2 = ODEProblem(eqn2, u02, tspan2, p2)
sol2 = solve(prob2)
plot(sol2, lab=["Position" "Velocity"])

# ---
# Straight Line Motion with Position-Dependent Force
f(x) = x + 1.0
function eqn3(du, u, p, t)
    m = p[1]
    du[1] = u[2]
    du[2] = f(u[1]) / m
end

mass3 = 10.0
p3 = [mass3]
u03 = [1.0, 0.0]
tspan3 = (0.0, 10.0)
prob3 = ODEProblem(eqn3, u03, tspan3, p3)
sol3 = solve(prob3)
plot(sol3, lab=["Position" "Velocity"])

