#!/usr/bin/julia -f

include("ode-common.jl")

immutable HarmonicPotential{T}
    omega::T
    center::T
end

call(h::HarmonicPotential, x) = h.omega^2 .* (x - h.center).^2

@inline _calc_d2_grid(n, d) = ifelse(n == 0, π^2 / 3,
                                     ifelse(n % 2 == 0,
                                            2, -2) / n^2) / (d^2)

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

grid_size = 162
grid_space = 0.02 * 150 / grid_size
x_omega = 5π

x_center = grid_size * grid_space / 2
psi_init = complex(exp(-linspace(-3 * x_center, 3 * x_center, grid_size).^2))

# println(psi_init)
# error()

h = HarmonicHamiltonian(x_omega, grid_space, x_center)

println("start")
@time t, y = solve_ode(0.0, psi_init, h, 0.2, 0.0001)

imshow(log(log(abs(y) + 1)))
colorbar()

println()
readline()
