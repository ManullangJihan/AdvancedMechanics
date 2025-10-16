# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using DifferentialEquations
using Plots
using Plots: plot, plot!

const g = 9.81
function eqn(du, u, p, t)
    m, drag = p
    x, y, vx, vy = u
    v = sqrt(vx^2 + vy^2)
    du[1] = vx
    du[2] = vy
    du[3] = -drag * x * v/m
    du[4] = -drag * y * v/m - g
end

function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[2]
end

function affect!(integrator)
    integrator.u[4] = -integrator.u[4]
end

function floor_aff!(int)
    int.p[2] = 0.0
	int.u[3] = 0.0
    int.u[4] = 0.0
end

## Set Parameters
# Ball's Mass 
m = 0.1

# Drag Coefficient
drag = 0.0002

# Time Final 
t_end = 2.0
is_bounced = true

## Set Initial Condition
vx0 = 5.0
vy0 = 3.0

p = [m, drag, g]
tspan = (0.0, t_end)
t = tspan[1]:0.01:t_end
x0 = 0.0
y0 = 0.0
u0 = [x0, y0, vx0, vy0]

if is_bounced
    global cb = ContinuousCallback(condition, affect!)
else
    global cb = ContinuousCallback(condition, floor_aff!)
end

prob = ODEProblem{true}(eqn, u0, tspan, p)
sol = Array(solve(prob, saveat=0.01, callback=cb))
ymax = maximum(sol[2, :]) + (maximum(sol[2, :]) * 0.2)

n_simulation = length(t)
anim = @animate for i in 1:n_simulation
    plot(sol[1, 1:i], sol[2, 1:i], xlabel="Distance (m)", y="Height (m)", label="Ball's Trajectory", xlims=(0.0, 10.0), ylims=(0.0, ymax), title="t=$(t[i])")
    scatter!([sol[1, i]], [sol[2, i]], label=false)
end

gif(anim)