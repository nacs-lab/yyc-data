#!/usr/bin/julia -f

# Quantum mechanical harmonic oscillator using pure positional space
# Fourier grid differentiation and RG ODE solver

include("ode-common.jl")

immutable HarmonicPotential{T}
    omega::T
    center::T
end

call(h::HarmonicPotential, x) = h.omega^2 .* (x - h.center).^2

@inline _calc_d2_grid(n, d) = ifelse(n == 0, π^2 / 3,
                                     ifelse(n % 2 == 0,
                                            2, -2) / n^2) / (d^2)

immutable Hamiltonian1D{T, P} <: ODEKernel
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

typealias HarmonicHamiltonian{To, Td} Hamiltonian1D{Td, HarmonicPotential{To}}

call(::Type{HarmonicHamiltonian}, omega, d, c) =
    Hamiltonian1D(d, HarmonicPotential(omega, c))

grid_size = 401
grid_space = 0.02 * 400 / (grid_size - 1)
x_omega = 5π

x_center = (grid_size + 1) * grid_space / 2
psi_init = complex(exp(-linspace(-2.5 * x_center, 1.5 * x_center, grid_size).^2))

h = HarmonicHamiltonian(x_omega, grid_space, x_center)

println("start")
@time t, y = solve_ode(0.0, psi_init, h, 0.2, 0.2 / 4000)

# 401 x 2000: stable, error -> 2.5e-7, 22s
# 401 x 4000: stable, error -> 0.8e-8, 42s

absy = abs(y)

figure()
imshow(absy[:, 1:2:end])
colorbar()

figure()
plot(absy[:, 1])
plot(absy[:, end])

figure()
plot(absy[:, 1] - absy[:, end])

println()
readline()
