# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using DifferentialEquations
using Plots

# Define the triangular mountain shape (x, y, height)
function mountain(x)
    return 0.2 * x.^2
end

# Define the ODE for the ball sliding down the mountain
function ball_sliding!(du, u, p, t)
    g = p[1]  # acceleration due to gravity
    
    # x and y coordinates
    x, y = u
    
    # Derivatives
    du[1] = sqrt(1 + (mountain'(x))^2)  # dx/dt = sqrt(1 + (dy/dx)^2)
    du[2] = du[1] * mountain'(x)  # dy/dt = (dy/dx) * dx/dt
    
    # Gravity component along the slope
    du[2] -= g * mountain'(x) / du[1]
end

# Initial conditions (starting position)
u0 = [0.0, mountain(0.0)]

# Parameters
g = 9.8  # acceleration due to gravity

# Set up the time span
tspan = (0.0, 5.0)

# Solve the ODE
prob = ODEProblem(ball_sliding!, u0, tspan, [g])
sol = solve(prob, Tsit5())

# Plot the mountain
x_values = range(-5, stop=5, length=100)
y_values = mountain.(x_values)

# Plot the ball sliding down the mountain
anim = @animate for i in 1:length(sol.t)
    plot(x_values, y_values, label="Mountain", xlabel="X", ylabel="Y", legend=:bottomright)
    plot!([sol.u[1][i]], [sol.u[2][i]], seriestype=:scatter, markersize=10, label="Ball")
    xlims!(-5, 5)
    ylims!(0, 2)
end

# To visualize the animation
gif(anim, "ball_sliding_animation.gif", fps=30)
