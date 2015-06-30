#!/usr/bin/julia -f

include("ode-common.jl")

immutable SimpleOscillator{T}
    omega::T
end

function call{T}(o::SimpleOscillator{T}, t, y, ydot)
    @inbounds ydot[1] = y[2]
    @inbounds ydot[2] = -o.omega^2 * y[1]
end
call{T}(o::SimpleOscillator{T}, t, y) = [y[2]; -o.omega^2 * y[1]]

@time t, y = solve_ode(0, [0.0, 1.0], SimpleOscillator(3.0), 10, 0.001)

plot(t, y'[:, 1])
plot(t, y'[:, 2])

println()
readline()
