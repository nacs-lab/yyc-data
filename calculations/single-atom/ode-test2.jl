#!/usr/bin/julia -f

include("ode-common.jl")

immutable HarmonicPotential{T}
    omega::T
    center::T
end

call(h::HarmonicPotential, x) = h.omega^2 .* (x - h.center).^2

@inline _calc_d2_grid(n, d) = ifelse(n == 0, π^2 / 3,
                                     ifelse(n % 2 == 0,
                                            2, -2) / n^2) / (d^2.5)

immutable Hamiltonian1D{T, P}
    d::T # grid spacing
    p::P # potential
end

function call(h::Hamiltonian1D, t, y, ydot)
    @inbounds for i in 1:length(y)
        x = i * h.d # coordinate
        v = h.p(x) * y[i] # potential term
        d = 0.0im # dynamic term
        for j in 1:length(y)
            d += _calc_d2_grid(i - j, h.d) * y[j]
        end
        ydot[i] = -im * (v + d)
    end
end

function call(h::Hamiltonian1D, t, y)
    ydot = similar(y)
    h(t, y, ydot)
    ydot
end

typealias HarmonicHamiltonian{To, Td} Hamiltonian1D{Td, HarmonicPotential{To}}

call(::Type{HarmonicHamiltonian}, omega, d, c) =
    Hamiltonian1D(d, HarmonicPotential(omega, c))

grid_size = 20
grid_space = 0.01
x_omega = 2π * 2

x_center = grid_size * grid_space / 2
psi_init = complex(exp(-linspace(-x_center, x_center, grid_size).^2))

h = HarmonicHamiltonian(x_omega, grid_space, x_center)

println("start")
# @time t, y = ODE.ode23(h, psi_init, linspace(0, 1, 400);
#                        abstol=3e-2, reltol=3e-2,
#                        initstep=0.02, points=:specified)
@time t, y = solve_ode(0.0, psi_init, h, 1, 0.001)

# @time ys = vector_hcat(y)

imshow(abs(y))

# plot(t, ys[:, 1])
# plot(t, ys[:, 2])

println()
readline()
