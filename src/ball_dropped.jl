# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using DifferentialEquations
using LinearAlgebra
using Plots
using Plots: plot, plot!

# ===
# Ball dropped

const g = 9.8
function eqn(du, u, p, t)
    du[1] = u[2]
    du[2] = -g
end

y0, v0 = 500, 0.0
u0 = [y0, v0]
tspan = (0.0, 10.0)
t = tspan[1]:0.1:tspan[2]

prob = ODEProblem(eqn, u0, tspan)
sol = solve(prob, saveat=t)
arr_sol = Array(sol)

ymin, ymax = minimum(arr_sol[1, :]), maximum(arr_sol[1, :])

anim = @animate for i in 1:length(t)
    scatter([1.0], [arr_sol[1, i]], xlabel="x", ylabel="Height", 
    xlims=(0.5, 1.5), ylims=(0, 550), markersize=10, label=false)
end
gif(anim, "gifs/ball_dropped.gif")

p = scatter(xlabel="x", ylabel="Height", xlims=(0.5, 1.5), ylims=(0, 550))
for i in 1:10:length(t)
    scatter!(p, [1.0], [arr_sol[1, i]], legend=false, color="blue")
end