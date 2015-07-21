#!/usr/bin/julia -f

# Quantum mechanical harmonic oscillator using split operator method with
# Float32 precision

immutable HarmonicPotential
    omega::Float32
    center::Float32
end

call(h::HarmonicPotential, x::Float32) = h.omega^2 .* (x - h.center).^2

immutable Hamiltonian1D{P}
    d::Float32 # grid spacing
    p::P # potential
end

function propagate(H::Hamiltonian1D, y0, t0::Float32, t1::Float32, dt::Float32)
    if t1 <= t0
        error("End time should be after start time.")
    end
    ts = t0:dt:t1
    nele = length(y0)
    nstep = length(ts)

    # Pre-allocate result and temporary arrays
    tmp = similar(y0)
    ys = Array{eltype(y0), 2}(nele, nstep)
    @inbounds ys[:, 1] = y0

    # FFT plan
    p_fft! = plan_fft!(tmp, flags=FFTW.MEASURE)

    # Propagators of x and p in it's diagonal form
    # The / 2 here is necessary to get intermediate results
    prop_x_2 = similar(y0)
    prop_p = similar(y0)

    k0 = Float32(2π / (nele * H.d))
    @inbounds for i in 1:nele
        x = i * H.d # coordinate
        prop_x_2[i] = exp(-im * H.p(x) * dt / 2)

        k1 = i - 1
        k2 = k1 - nele
        k = ifelse(k1 + k2 <= 0, k1, k2) * k0
        prop_p[i] = exp(-im * k^2 * dt)
    end

    @inbounds for i in 2:nstep
        for j in 1:nele
            tmp[j] = ys[j, i - 1] * prop_x_2[j]
        end
        p_fft! * tmp
        for j in 1:nele
            tmp[j] *= prop_p[j]
        end
        p_fft! \ tmp
        for j in 1:nele
            ys[j, i] = tmp[j] * prop_x_2[j]
        end
    end
    collect(ts), ys
end

typealias HarmonicHamiltonian Hamiltonian1D{HarmonicPotential}

call(::Type{HarmonicHamiltonian}, omega, d, c) =
    Hamiltonian1D(d, HarmonicPotential(omega, c))

grid_size = 2048
grid_space = 0.01f0 # * 1000 / (grid_size - 1)
x_omega = 1f0π

x_center = (grid_size + 1) * grid_space / 2
x_center2 = x_center
psi_init = complex(exp(-linspace(-2.0f0 * x_center2,
                                 1.0f0 * x_center2, grid_size).^2))

h = HarmonicHamiltonian(x_omega, grid_space, x_center)

println("start")
@time t, y = propagate(h, psi_init, 0.0f0, 0.2f0, 0.2f0 / 100)
gc()
@time t, y = propagate(h, psi_init, 0.0f0, 100.0f0, 1.0f0 / 100)

# exit()

using PyPlot

absy = abs(y)

figure()
# imshow(absy[:, 1:2:end])
imshow(absy)
colorbar()

figure()
imshow(log(log1p(absy[:, 10:end])))
colorbar()

figure()
plot(absy[:, 1])
plot(absy[:, end])

figure()
plot(absy[:, 1] - absy[:, end])

show()
