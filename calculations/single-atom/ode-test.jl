#!/usr/bin/julia -f

include("ode-common.jl")

immutable SimpleOscillator{T}
    omega::T
end

call{T, V}(o::SimpleOscillator{T}, t, y::V) = [y[2]; -o.omega^2 * y[1]]

t, y = ode23s(SimpleOscillator(3.0), [0.0; 1.0], linspace(0, 10, 10000))

ys = vector_hcat(y)

plot(t, ys[:, 1])
plot(t, ys[:, 2])

println()
readline()
