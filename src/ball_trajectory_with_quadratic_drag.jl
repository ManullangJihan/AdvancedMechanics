using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

Cd = 0.4
ro = 1.225
ball_radius = 1.0
m = 1.0

A = π * ball_radius^2
K = 0.5*Cd*ro*A
g = 9.81

function eqn(du, u, p, t)

    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    du[4] = (-K/m) * u[4] * sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    du[5] = (-K/m) * u[5] * sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    du[6] = -g - (K/m) * u[6] * sqrt(u[4]^2 + u[5]^2 + u[6]^2)
end

tspan = (0.0, 10.0)
u0 = [0.0, 0.0, 0.0, 2.0, 2.0, 2.0]
prob = ODEProblem(eqn, u0, tspan)
sol = solve(prob)
arr_sol = Array(sol)
scatter(arr_sol[1, :], arr_sol[2, :], arr_sol[3, :])