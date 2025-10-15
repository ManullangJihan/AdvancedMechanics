# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")


using DifferentialEquations
using Plots
using Plots: plot, plot!


function suspension_sys!(du, u, p, t)
    m1, m2, k1, k2, b1, b2 = p
    du[1] = u[2]
    du[2] = -(b1*b2)/(m1*m2)*u[1] + ((b1/m1)*((b1/m1)+(b1/m2)+(b2/m2)))-(k1/m1)*u[3] - (b1/m1)*u[4]
    du[3] = b2/m2*u[1] - ((b1/m1)+(b1/m2)+(b2/m2))*u[3] + u[4]
    du[4] = k2/m2*u[1] - ((k1/m1)+(k1/m2)+(k2/m2))*u[3]
end


m1 = 2500
m2 = 320
k1 = 80000
k2 = 500000
b1 = 350
b2 = 15020

p = (m1, m2, k1, k2, b1, b2)

x1, v1, x2, v2 = 4.0, 0.0, 2.0, 0.0
u0 = [x1, x2, v1, v2]
t_end = 10.0
tspan = (0.0, t_end)
t = tspan[1]:0.1:tspan[2]
prob = ODEProblem(suspension_sys!, u0, tspan, p)
sol = solve(prob, saveat=t)
arr_sol = Array(sol)

x1s = arr_sol[1, :]
x2s = arr_sol[2, :]

plot(t, arr_sol[1, :])
plot(t, arr_sol[3, :])

l1 = 5.0
l2 = 5.0
W1 = 2.0
W2 = 2.0
H1 = 3.0
H2 = 1.0
x10 = H2 + l2 + l1

rectangle(x, y, w, h) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

x1_ = x1s[1]
x2_ = x2s[1]
p = plot([W2/2, W2/2], [0.0, x2_])
plot!(p, rectangle(0.0, x2_, W2, H2))

plot!(p, [W2/2, W2/2], [x2_ + H2, x1_])
plot!(p, rectangle(0.0, x1_, W1, H1))


p = plot(rectangle(0.0, x1s[1] + 1.0, W1, H1), legend=false, xlims=(0, 20), ylims=(0.0, 20.0))
plot!(p, [2.5, 2.5],  [l2 + H2 + x1s[1], x1s[1]], linewidth=3.0)

plot!(p, rectangle(W1/2 - 1.0, l2 + x2s[1], W2, H2), xlims=(-l1, 15.0), ylims=(0.0, 6.0), legend=false)
plot!(p, [2.5, 2.5], [l1, x2s[1]], linewidth=3.0)