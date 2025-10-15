using LinearAlgebra
using ControlSystems
using Plots
using Plots: plot, plot!
using DifferentialEquations

m = 1
M = 5
L = 2
g = -9.8
d = 1

b = 1 # pendulum down (b=-1)

A = [0 1 0 0;
     0 -d/M b*m*g/M 0;
     0 0 0 1;
     0 -b*d/(M*L) -b*(m+M)*g/(M*L) 0]
B = [0; 1/M; 0; b*1/(M*L)]
C = [1 0 0 0]
D = 0
sys = ss(A, B, C, D)
Q  = [
     1 0 0 0; 
     0 1 0 0; 
     0 0 1 0; 
     0 0 0 1]    # Weighting matrix for state
R = 0.1
K = lqr(sys, Q, R)
ref = [1, 0, π, 0];
u(x, t)  = -K*(x .- wr)    

# Simulation
u(x, t)  = -K * (x .- ref);
function pendcart(dx, x, p, t)
          Sx = sin(x[3])
          Cx = cos(x[3])
          D = m*L*L*(M+m*(1-Cx^2))
          
          dx[1] = x[2]
          dx[2] = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x[4]^2*Sx - d*x[2])) + m*L*L*(1/D)*u(x, t)[1]
          dx[3] = x[4]
          dx[4] = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x[4]^2*Sx - d*x[2])) - m*L*Cx*(1/D)*u(x, t)[1]
end

tspan = (0.0, 10.0)
t = tspan[1]:0.1:tspan[2]
x0 = [3, 0, π-0.2, 0]
prob = ODEProblem(pendcart, x0, tspan)
sol = solve(prob, Tsit5(), saveat=t)
x_res = Array(sol)

# dimensions
# dimensions
W = 1 * sqrt(M/5)   # cart width
H = 0.5 * sqrt(M/5) # cart height
wr = 0.1            # wheel radius
mr = 0.3*sqrt(m)    # mass radius

# positions
y = wr/2 + H/2      # cart vertical position

rectangle(x, y, w, h) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
function circleShape(h, k, r)
     θ = LinRange(0, 2π, 500)
     h .+ r*sin.(θ), k .+ r*cos.(θ)
end

function drawpend(state, i)
     x = state[1]
     θ = state[3]
     pendx = x + L*sin(θ)
     pendy = y - L*cos(θ)

     plot(rectangle(x-W/2, y-H/2, W, H), xlim=(0.0, 5), ylim=(-2, 2.5), size=(1200, 800))
     plot!(circleShape(x - 0.9 * W/2, 0, wr), seriestype=[:shape], lw=0.5, c=:blue, 
          linecolor=:black, legend=false, fillalpha=0.2, aspect_ratio=1
     )
     plot!(circleShape(x + 0.9 * W/2, 0, wr), seriestype=[:shape], lw=0.5, c=:blue, 
          linecolor=:black, legend=false, fillalpha=0.2, aspect_ratio=1
     )
     plot!([x, pendx], [y, pendy], color=:black, linewidth=3)

     plot!(circleShape(pendx, pendy, 0.2), seriestype=[:shape], color=:yellow, 
     linecolor=:black, fillalpha=1.0, legend=false)
end

t_simulation = size(x_res, 2)
anim = @animate for i in 1:t_simulation
     state = x_res[:, i]
     drawpend(state, i)
end

gif(anim, "inverted_pendulum.gif")