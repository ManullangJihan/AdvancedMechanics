


# ---
# Particle Sliding on a Parabola
function eqn8(du, u, p, t)
    μ = p[1]
    du[1] = u[2]
    du[2] = -2*μ*u[1] * (g+2*μ*u[2]^2) / (1+4*μ^2 * u[1]^2)
end

m8 = 1.0
μ8 = 1.0
p8 = [μ8]
u08 = [1.0, 0.0]
tspan8 = (0.0, 10.0)
t8 = tspan8[1]:0.1:tspan8[2]
prob8 = ODEProblem(eqn8, u08, tspan8, p8)
sol8 = solve(prob8, saveat=t8)
arr_sol8 = Array(sol8)
x8, ẋ8 = arr_sol8[1, :], arr_sol8[2, :]
N = @. m8 * (g + 2 * μ8 * ẋ8^2) / sqrt(1 + 4 * μ8^2 * x8^2)
plot(sol8, lab=["Position" "Velocity"])

# ---
# A Particle Sliding on a Hemisphere

function eqn9(du, u, p, t)
    R = p[1]
    du[1] = u[2]
    du[2] = -g/R * sin(u[1])
end

R9 = 10.0
p9 = [R9]
u09 = [0.0, 1.0]
tspan9 = (0.0, 10.0)
t9 = tspan9[1]:0.1:tspan9[2]
prob9 = ODEProblem(eqn9, u09, tspan9, p9)
sol9 = solve(prob9, saveat=t9)
arr_sol9 = Array(sol9)
z9, zdot9 = arr_sol9[1, :], arr_sol9[2, :]
theta9 = z9 .- 90 
x9 = R9 .* cos.(theta9)
y9 = R9 .* sin.(theta9)

plot(x9, y9)

# ---
# Particle Sliding on Parabola
using DifferentialEquations
using Plots
using Plots: plot, plot!

g = 9.8
function eqn(du, u, p, t)
    μ = p[1]
    du[1] = u[2]
    du[2] = -2μ * u[1] * (g + 2μ * u[2]^2) / (1 + 4μ^2 * u[1]^2)
end

μ = 1.0
p = [μ]
u0 = [1.0, 0.0]
t_end = 10.0
tspan = (0.0, t_end)
t = tspan[1]:0.1:tspan[2]

prob = ODEProblem(eqn, u0, tspan, p)
sol = solve(prob, saveat=t)
arr_sol = Array(sol)
y(x) = μ * x^2
position = arr_sol[1, :]
xs = -1.5:0.1:1.5

anim = @animate for x in position
    p = plot(xs, y.(xs), legend=false)
    scatter!(p, [x], [y.(x)], label="")
end

gif(anim, "particle_sliding.gif")