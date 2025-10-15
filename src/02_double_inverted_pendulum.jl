# https://chat.openai.com/share/3d2e0c3f-b9ac-42e3-a798-e488095e85d0

using DifferentialEquations
using Plots
using Plots: plot, plot!

const g = 9.81
function double_pendulum!(du, u, p, t)
    # Parameters
    m1, m2, l1, l2, c1, c2 = p

    # Extract variables
    θ1, θ2, ω1, ω2 = u

    # Equations of motion
    du[1] = ω1
    du[2] = ω2

    Δθ = θ2 - θ1
    den1 = (m1 + m2) * l1 - m2 * l1 * cos(Δθ)^2
    den2 = (m1 + m2) * l2 - m2 * l2 * cos(Δθ)^2

    du[3] = (m2 * l2 * ω2^2 * sin(Δθ) * cos(Δθ) +
             m2 * g * sin(θ2) * cos(Δθ) +
             m2 * l2 * ω2^2 * sin(Δθ) -
             (m1 + m2) * g * sin(θ1) - c1 * ω1) / den1

    du[4] = ((m1 + m2) * l1 * ω1^2 * sin(Δθ) * cos(Δθ) -
             (m1 + m2) * g * sin(θ2) * cos(Δθ) +
             (m1 + m2) * g * sin(θ1) * cos(Δθ) +
             m2 * l2 * ω2^2 * sin(Δθ) - c2 * ω2) / den2
end

# parameters
m1, m2 = 1.0, 1.0
l1, l2 = 1.0, 1.0

# Damping coefficients
c1, c2 = 0.1, 0.1

p = [m1, m2, l1, l2, c1, c2]

# Initial conditions
u0 = [π/4, π/4, 0.0, 0.0]

# Time span
tspan = (0.0, 10.0)
t = tspan[1]:0.1:tspan[2]

# Solve the system of ODEs
prob = ODEProblem(double_pendulum!, u0, tspan, p)
sol = solve(prob, Vern7(), saveat=t)
arr_sol = Array(sol)

θ1, θ2 = arr_sol[1, :], arr_sol[2, :]

x1 = l1 .* sin.(θ1)
y1 = -l1 .* cos.(θ1)
x2 = l1 .* sin.(θ1) .+ l2 .* sin.(θ2)
y2 = -l1 .* cos.(θ1) .- l2 .* cos.(θ2)

function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

anim = @animate for (x1_, y1_, x2_, y2_) in zip(x1, y1, x2, y2)

    plot([0.0, x1_], [0.0, y1_], xlim=(-2.0, 2.0), ylim=(-3.5, 1.0), linewidth=5.0, xaxis=false, yaxis=false)
    plot!(
        circleShape(x1_, y1_, 0.2), seriestype=[:shape,], lw=0.5, c=:blue, 
        linecolor=:black, legend=false, fillalpha=0.2, aspect_ratio=1
    )
    hline!([0.0], linewidth=3.0)

    plot!([x1_, x2_], [y1_, y2_],linewidth=5.0, xaxis=false, yaxis=false)
    plot!(
        circleShape(x2_, y2_, 0.2), seriestype=[:shape,], lw=0.5, c=:blue, 
        linecolor=:black, legend=false, fillalpha=0.2, aspect_ratio=1
    )
end

gif(anim, "double_pendulum.gif")
