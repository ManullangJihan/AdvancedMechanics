
const G = 6.674 * 10e−11
const Re = 6.371 * 10^6

function eqn(du, u, p, t)
    m, c = p
    y, v = u
    du[1] = u[2]
    du[2] = - g + c/m * v^2
end

m = 10.0
c = 0.05
p = [m, c]
u0 = [1000.0, 0.0]
tspan = (0.0, 20.0)

prob = ODEProblem(eqn, u0, tspan, p)
sol = solve(prob)
plot(sol, lab=["Position" "Velocity"])