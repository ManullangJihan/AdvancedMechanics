# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

const g = 9.81
function airdrag(du, u, p, t)
    m, c, w = p
    du[1] = u[2]
    du[2] = -c/m * u[2] * sqrt(u[2]^2 + u[4]^2)
    du[3] = u[4]
    du[4] = -w/m - c/m * u[4] * sqrt(u[2]^2 + u[4]^2)
end

m = 10
c = 0.2
w = m * g
p = [m, c, w]
tspan = (0.5, 5.0)
u₀ = [0.0, 100., 0., 10.]

function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[3]
end

function affect!(integrator)
    integrator.u[4] = -integrator.u[4]
end

cb = ContinuousCallback(condition, affect!)

problem = ODEProblem(airdrag, u₀, tspan, p)
sol = solve(problem, saveat=0.1, callback = cb)
t = round.(sol.t, digits=2)
sol = Array(sol)
ymin, ymax = minimum(sol[3, :]), maximum(sol[3, :]) + maximum(sol[3, :])*0.1
xmin, xmax = minimum(sol[1, :]), maximum(sol[1, :]) + maximum(sol[1, :])*0.1

anim = @animate for i in 1:length(t)
    plot(sol[1, 1:i], sol[3, 1:i], xlabel="Distance x", ylabel="Height",
    title="Ball Trajectory at time: $(t[i])", xlims=(xmin, xmax), ylims=(ymin, ymax),
    legend=false, linewidth=5)
    scatter!([sol[1, i]], [sol[3, i]])
end

gif(anim)