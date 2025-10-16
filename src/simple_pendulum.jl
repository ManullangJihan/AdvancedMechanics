# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

const b = 0.1
const g = -9.8
const L = 1.0
const m = 2.0

function eqn(du, u, p, t)
    du[1] = u[2]
    du[2] = -g/L * sin(u[1]) - b/m * u[2]
end

x1 = 1.1
x2 = 0.0
u0 = [x1, x2]
tspan = (0.0, 10.0)
t = tspan[1]:0.1:tspan[2]

prob = ODEProblem(eqn, u0, tspan)
sol = solve(prob, saveat=t)
arr_sol = Array(sol)
Plots.plot(sol)

xs = sin.(arr_sol[1, :]) * L
ys = cos.(arr_sol[1, :]) * L

function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

anim = @animate for (x, y) in zip(xs, ys)
    plot([0.0, x], [0.0, y], xlim=(-2.0, 2.0), ylim=(-1.5, 1.0), linewidth=5.0, xaxis=false, yaxis=false)
    plot!(circleShape(x, y, 0.2), seriestype=[:shape,], lw=0.5, c=:blue, 
    linecolor=:black, legend=false, fillalpha=0.2, aspect_ratio=1)
    hline!([0.0], linewidth=3.0)
end

gif(anim, "simple_pendulum.gif")
