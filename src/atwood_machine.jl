# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using DifferentialEquations
using Plots
using Plots: plot, plot!

const g = 9.8
function eqn(du, u, p, t)
    x, v = u
    m1, m2 = p
    du[1] = dx = v
    du[2] = dv = (m2 - m1)*g/(m2+m1) 
end

x₀, v₀ = 0.0, 0.0
u₀ = [x₀, v₀]
m1, m2 = 1.1, 1.0
p = (m1, m2)

tspan = (0.0, 1.5)
prob = ODEProblem(eqn, u₀, tspan, p)
sol = solve(prob, saveat=0.1)
t = sol.t
sol = Array(sol)
N = length(t)

function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

x = sol[1, :]

anim = @animate for i in 1:N
    # plot([-1.0, -1.0], [0.0, -1*(-x[i]+1.2)], legend=false, 
    # xticks=false, yticks=false,
    # xlims=(-1.2, 0.5), ylims=(-2.5, 0.2), linewidth=3, color=:blue)
    plot([-1.0, -1.0], [0.0, -1*(-x[i]+1.2)], legend=false, 
    xticks=false, yticks=false,
    xlims=(-1.5, 1.5), ylims=(-2.5, 0.2), linewidth=3, color=:blue)
    plot!(circleShape(-1.0, -1*(-x[i]+1.2), 0.09), seriestype=[:shape])
    plot!([-1.0, 0.0], [0.0, 0.0], linewidth=3, color=:blue)

    plot!([0.0, 0.0], [0.0, -1*(x[i]+1)], linewidth=3, color=:blue)
    plot!(circleShape(0.0, -1*(x[i]+1), 0.09), seriestype=[:shape])
end
gif(anim, "gifs/atwood_machine.gif")