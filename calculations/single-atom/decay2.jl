#!/usr/bin/julia -f

# Decay of a driven two level system with Monte Carlo

using PyPlot

immutable StationaryAtom
    Γ::Float64
    Ω::Float64
    δ::Float64
end

immutable StationaryPropagator
    H::StationaryAtom
    dt::Float64
end

@inline function propagate(P::StationaryPropagator, ψ0, ψs::Matrix{Complex128})
    nstep = size(ψs, 2)
    @inbounds ψs[:, 1] = ψ0

    H = P.H
    dt = P.dt

    next_decay = -log(rand()) / H.Γ

    cos_dt = cos(H.Ω * dt)
    sin_dt = sin(H.Ω * dt)

    @inbounds for i in 2:nstep
        t = i * P.dt
        if t > next_decay
            next_decay += -log(rand()) / P.H.Γ
            ψ_e = Complex128(0)
            ψ_g = Complex128(1)
        else
            ψ_e = ψs[2, i - 1]
            ψ_g = ψs[1, i - 1]
        end

        # Hamiltonian is
        # H = [0, Ω exp(im δ t)
        #      Ω exp(-im δ t), 0]
        #   = Ω * (cos(δ t) σ_x + sin(δ t) σ_y)

        θ = H.δ * (t + dt / 2)
        cos_t = cos(θ)
        sin_t = sin(θ)

        # Propagator is
        # exp(im H Δt) = exp(im * Ω * (cos_t σ_x + sin_t σ_y) * Δt)
        #              = cos(Ω Δt) + im * (cos_t σ_x + sin_t σ_y) * sin(Ω Δt)
        #              = cos_dt + im cos_t sin_dt σ_x + im sin_t sin_dt σ_y
        #              = [cos_dt, im (cos_t + im sin_t) sin_dt
        #                 im (cos_t - im sin_t) sin_dt, cos_dt]

        T22 = T11 = cos_dt
        T12 = (im * cos_t - sin_t) * sin_dt
        T21 = (im * cos_t + sin_t) * sin_dt

        ψs[2, i] = T11 * ψ_e + T12 * ψ_g
        ψs[1, i] = T22 * ψ_g + T21 * ψ_e
    end
    ψs
end

@inline function propagate(H::StationaryAtom, ψ0, dt, ψs::Matrix{Complex128})
    P = StationaryPropagator(H, dt)
    propagate(P, ψ0, ψs)
end

@inline function propagate(H::StationaryAtom, ψ0, dt, nstep::Integer)
    # Pre-allocate result and temporary arrays
    ψs = Array{Complex{Float64}, 2}(2, nstep)
    propagate(H, ψ0, dt, ψs)
end

function propagate_average(H, dt, nstep, ntime)
    ψ_init = Complex{Float64}[0, 1]
    ψs = Matrix{Complex{Float64}}(2, nstep)
    ps = zeros(Float64, (2, nstep))
    P = StationaryPropagator(H, dt)
    @inbounds for i in 1:ntime
        propagate(P, ψ_init, ψs)
        for j in 1:nstep
            ps[1, j] += abs2(ψs[1, j])
            ps[2, j] += abs2(ψs[2, j])
        end
    end
    @inbounds for j in 1:nstep
        ps[1, j] /= ntime
        ps[2, j] /= ntime
    end
    ps
end

# Γ, Ω, δ
h = StationaryAtom(10e6, 40e6, 40e6)

println("start")
@time ps = propagate_average(h, 1e-9, 10, 1)
gc()
@time ps = propagate_average(h, 1e-9, 500, 100_000)

figure()
plot(ps'[:, 1])
plot(ps'[:, 2])

show()
