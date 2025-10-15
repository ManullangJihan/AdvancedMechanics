# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using Symbolics
using DifferentialEquations
using Plots
using Plots: plot, plot!

@variables t, g, m1, m2, l1, l2
@variables θ₁, θ₂, ω₁, ω₂
D = Differential(t)
x1 = l1 * sin(θ₁)
y1 = -l1 * cos(θ₁)
x2 = l1 * sin(θ₁) - l2 * sin(θ₂)
y2 = -l1 * cos(θ₁) - l2 * cos(θ₂)

# Kinetic
theta1_d = D(θ₁)
theta2_d = D(θ₂)
theta1_dd = D(theta1_d)
theta2_dd = D(theta2_d)

