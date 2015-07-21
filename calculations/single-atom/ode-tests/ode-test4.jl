#!/usr/bin/julia -f

# Quantum mechanical harmonic oscillator using FFT differentiation and
# RG ODE solver

include("ode-common.jl")

immutable HarmonicPotential{T}
    omega::T
    center::T
end

call(h::HarmonicPotential, x) = h.omega^2 .* (x - h.center).^2

immutable Hamiltonian1D{T, P} <: ODEKernel
    d::T # grid spacing
    p::P # potential
end

function call(h::Hamiltonian1D, t, y, ydot)
    len = length(y)
    @inbounds for i in 1:len
        ydot[i] = y[i]
    end
    y_fft! * ydot
    k0 = 2π / (len * h.d)
    @inbounds for i in 1:len
        k1 = i - 1
        k2 = k1 - len
        k = ifelse(k1 + k2 <= 0, k1, k2) * k0
        ydot[i] = ydot[i] * k^2
    end
    y_ifft! * ydot
    @inbounds for i in 1:len
        x = i * h.d # coordinate
        v = h.p(x) * y[i] # potential term
        ydot[i] = -im * (v + ydot[i])
    end
end

typealias HarmonicHamiltonian{To, Td} Hamiltonian1D{Td, HarmonicPotential{To}}

call(::Type{HarmonicHamiltonian}, omega, d, c) =
    Hamiltonian1D(d, HarmonicPotential(omega, c))

grid_size = 1001
grid_space = 0.02 * 1000 / (grid_size - 1)
x_omega = 5π

x_center = (grid_size + 1) * grid_space / 2
psi_init = complex(exp(-linspace(-2.5 * x_center, 1.5 * x_center, grid_size).^2))

const y_fft! = plan_fft!(copy(psi_init), flags=FFTW.MEASURE)
const y_ifft! = plan_ifft!(copy(psi_init), flags=FFTW.MEASURE)

h = HarmonicHamiltonian(x_omega, grid_space, x_center)

println("start")
@time t, y = solve_ode(0.0, psi_init, h, 1.0, 0.2 / 4000)
gc()
@time t, y = solve_ode(0.0, psi_init, h, 1.0, 0.2 / 4000)

# exit()

absy = abs(y)

figure()
imshow(absy[:, 1:2:end])
colorbar()

figure()
plot(absy[:, 1])
plot(absy[:, end])

figure()
plot(absy[:, 1] - absy[:, end])

show()
