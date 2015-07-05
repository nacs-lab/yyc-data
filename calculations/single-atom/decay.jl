#!/usr/bin/julia -f

# Decay of a driven two level system with Monte Carlo

using PyPlot

immutable StationaryAtom
    Γ::Float64
    Ω::Float64
    δ::Float64
end

function propagate(H::StationaryAtom, ψ0, dt, ψs::Matrix{Complex128})
    @assert length(ψ0) == 2

    nstep = size(ψs, 2)
    @inbounds ψs[:, 1] = ψ0

    # Hamiltonian is
    # H = [δ / 2, Ω
    #      Ω, -δ / 2] = (δ / 2) σ_z + Ω σ_x = c_σ * (c_σz σ_z + c_σx σ_x)

    c_σ = √(H.δ^2 / 4 + H.Ω^2)
    c_σx = H.Ω / c_σ
    c_σz = H.δ / 2 / c_σ

    # Propagator is
    # exp(im H Δt) = exp(im * c_σ * (c_σz σ_z + c_σx σ_x) * Δt)
    #              = cos(c_σ Δt) + im * (c_σz σ_z + c_σx σ_x) * sin(c_σ Δt)
    #              = cos_dt + (im c_σz sin_dt) σ_z + (im c_σx sin_dt) σ_x
    #              = [cos_dt + im c_σz sin_dt, im c_σx sin_dt
    #                 im c_σx sin_dt         , cos_dt - im c_σz sin_dt]

    cos_dt = cos(c_σ * dt)
    sin_dt = sin(c_σ * dt)

    T11 = cos_dt + im * c_σz * sin_dt
    T22 = cos_dt - im * c_σz * sin_dt
    T21 = T12 = im * c_σx * sin_dt

    next_decay = -log(rand()) / H.Γ

    @inbounds for i in 2:nstep
        t = i * dt
        if t > next_decay
            next_decay += -log(rand()) / H.Γ
            ψ_e = Complex128(0)
            ψ_g = Complex128(1)
        else
            ψ_e = ψs[2, i - 1]
            ψ_g = ψs[1, i - 1]
        end

        ψs[2, i] = T11 * ψ_e + T12 * ψ_g
        ψs[1, i] = T22 * ψ_g + T21 * ψ_e
    end
    ψs
end

function propagate(H::StationaryAtom, ψ0, dt, nstep::Integer)
    # Pre-allocate result and temporary arrays
    ψs = Array{Complex{Float64}, 2}(2, nstep)
    propagate(H, ψ0, dt, ψs)
end

function propagate_average(H, dt, nstep, ntime)
    ψ_init = Complex{Float64}[0, 1]
    ψs = Matrix{Complex{Float64}}(2, nstep)
    ps = zeros(Float64, (2, nstep))
    @inbounds for i in 1:ntime
        propagate(H, ψ_init, dt, ψs)
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
h = StationaryAtom(20e6, 40e6, 20e6)

println("start")
@time ps = propagate_average(h, 1e-9, 10, 1)
gc()
@time ps = propagate_average(h, 5e-10, 500, 100_000)

figure()
plot(ps'[:, 1])
plot(ps'[:, 2])

println()
readline()
