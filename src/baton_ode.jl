# https://www.mathworks.com/help/matlab/math/solve-equations-of-motion-for-baton.html

using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

m1 = 0.1
m2 = 0.1
L = 2.0
g = 9.81

function eqn(du, u, p, t)
    x, xv, y, yv, θ, θv = u
    du[1] = dx = xv
    du[2] = dxv = m2*L*θv^2*cos(θ)
    du[3] = dy = yv
    du[4] = dyv = m2*L*θv^2*sin(θ)-(m1+m2)*g
    du[5] = dθ = θv
    du[6] = dθv = -g*L*cos(θ)
end

tspan = (0.0, 20.0)
u0 = [0, 4.0, L, 25.0, -π, 2]
problem = ODEProblem(eqn, u0, tspan)
sol = solve(problem)
t = sol.t
arr_sol = Array(sol)

# TODO: Make animation

anim = @animate for j = 1:length(t)
    theta = arr_sol[5, j]
    X = arr_sol[1, j]
    Y = arr_sol[3, j]
    xvals = [X, X+L*cos(theta)]
    yvals = [Y, Y+L*sin(theta)]
    plot(xvals, yvals, xlims=(0.0, 50.0), ylims=(0.0, 80.0), legend=false, lw=3,
    xaxis=false, yaxis=false)
    scatter!([xvals[1]], [yvals[1]], color=:red)
    scatter!([xvals[2]], [yvals[2]], color=:green)
end

gif(anim, "gifs/baton_system.gif")
