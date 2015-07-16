#!/usr/bin/julia -f

# 1D harmonic trap (non-magic wavelength) using Monte Carlo

##
# Hamiltonian and jump operators

immutable HTrap{T}
    m::T
    ω::NTuple{2, T}
end

@inline function get_potential(h::HTrap, x, idx)
    ω = h.ω[idx]
    h.m * ω^2 * x^2 / 2
end
@inline function get_kinetic(h::HTrap, k)
    k^2 / (2 * h.m)
end

immutable OpticalDecay{T}
    k::T
    Γ::T
end

immutable OpticalDrive{T}
    k::T
    Ω::T
    δ::T
    τ_θ::T
end

immutable HSystem{H, T, N}
    trap::H
    decay::OpticalDecay{T}
    drives::NTuple{N, OpticalDrive{T}}
end

type PhaseTracker{T}
    drive::OpticalDrive{T}
    phase::T
    prev_t::T

    total_phase::T
    sin_t::T
    cos_t::T
end

function call{T}(::Type{PhaseTracker}, drive::OpticalDrive{T})
    PhaseTracker{T}(drive, T(0), T(0), T(0), T(0), T(0))
end

function phase_tracker_init{T}(track::PhaseTracker{T})
    if isfinite(track.drive.τ_θ)
        track.phase = rand(T) * 2π
    else
        track.phase = 0
    end
    track.prev_t = 0
end

function phase_tracker_next{T}(track::PhaseTracker{T}, t::T)
    prev_t = track.prev_t
    if prev_t >= t
        return track.phase
    end
    track.prev_t = t
    if !isfinite(track.drive.τ_θ)
        return track.phase
    end
    δt = (t - prev_t) / track.drive.τ_θ
    δθ = sqrt(δt) * (rand(T) - 0.5) * π
    track.phase = (track.phase + δθ) % 2π
    track.phase
end

function phase_tracker_update{T}(track::PhaseTracker{T}, t::T)
    phase = phase_tracker_next(track, t)
    track.total_phase = (phase + track.drive.δ * t) % 2π
    track.sin_t = sin(track.total_phase)
    track.cos_t = cos(track.total_phase)
    nothing
end

##
# Propagator

immutable SystemPropagator{H, T, N, P, PI}
    H::HSystem{H, T, N}
    dt::T
    dx::T
    nstep::Int
    nele::Int

    tmp::Matrix{Complex{T}} # 2 x nele
    p_fft!::P
    p_bfft!::PI

    # Energy
    E_k::Vector{T} # nele
    E_x::NTuple{2, Vector{T}} # nele

    # Propagator
    P_k::Vector{Complex{T}} # nele
    P_x2::NTuple{2, Vector{Complex{T}}} # nele

    # Phase factors
    sin_decay::Vector{T} # nele
    cos_decay::Vector{T} # nele

    sin_drive::NTuple{N, Vector{T}} # nele
    cos_drive::NTuple{N, Vector{T}} # nele

    drive_phase::NTuple{N, PhaseTracker{T}}
end

function get_k_sin_cos{T}(::Type{T}, k, dx, nele)
    sins = Vector{T}(nele)
    coss = Vector{T}(nele)
    nele_2 = nele ÷ 2
    @inbounds for i in 1:nele
        x = (i - nele_2) * dx # coordinate
        θ = k * x
        sins[i] = sin(θ)
        coss[i] = cos(θ)
    end
    sins, coss
end

function get_drive_sin_cos{H, T}(h::HSystem{H, T}, dx, nele, idx)
    drive = h.drives[idx]
    k = drive.k
    get_k_sin_cos(T, k, dx, nele)
end

@generated function get_drives_sin_cos{H, T, N}(h::HSystem{H, T, N}, dx, nele)
    body = Expr(:block)
    sins = Expr(:tuple)
    coss = Expr(:tuple)
    resize!(body.args, N + 1)
    resize!(sins.args, N)
    resize!(coss.args, N)
    @inbounds for i in 1:N
        @gensym sym_sin
        @gensym sym_cos
        expr = quote
            $sym_sin, $sym_cos = get_drive_sin_cos(h, dx, nele, $i)
        end
        body.args[i] = expr
        sins.args[i] = sym_sin
        coss.args[i] = sym_cos
    end
    body.args[N + 1] = Expr(:tuple, sins, coss)
    body
end

@generated function get_drives_tracker{H, T, N}(h::HSystem{H, T, N})
    trackers = Expr(:tuple)
    resize!(trackers.args, N)
    @inbounds for i in 1:N
        trackers.args[i] = :(PhaseTracker(h.drives[$i]))
    end
    trackers
end

function SystemPropagator{H, T, N}(h::HSystem{H, T, N}, dt::T, dx::T,
                                   nstep, nele)
    tmp = Matrix{Complex{T}}(2, nele) # 2 x nele

    # FFT plan
    p_fft! = plan_fft!(tmp, 2, FFTW.MEASURE)
    p_bfft! = plan_bfft!(tmp, 2, FFTW.MEASURE)
    # p_fft! = plan_fft!(tmp, 2, flags=FFTW.MEASURE)
    # p_bfft! = plan_bfft!(tmp, 2, flags=FFTW.MEASURE)

    E_k = Vector{T}(nele)
    E_xg = Vector{T}(nele)
    E_xe = Vector{T}(nele)

    P_k = Vector{Complex{T}}(nele)
    P_x2g = Vector{Complex{T}}(nele)
    P_x2e = Vector{Complex{T}}(nele)

    k0 = 2π / (nele * dx)
    nele_2 = nele ÷ 2
    @inbounds for i in 1:nele
        x = (i - nele_2) * dx
        e_xg = get_potential(h.trap, x, 1)
        e_xe = get_potential(h.trap, x, 2)
        E_xg[i] = e_xg
        E_xe[i] = e_xe
        P_x2g[i] = exp(im * e_xg * dt / 2)
        P_x2e[i] = exp(im * e_xe * dt / 2)

        k1 = i - 1
        k2 = k1 - nele
        k = ifelse(k1 + k2 <= 0, k1, k2) * k0

        e_k = get_kinetic(h.trap, k)
        E_k[i] = e_k
        P_k[i] = exp(im * e_k * dt)
    end

    # Precompute the spacial phase factors of the light
    # The phase factor from time can be calculated at each time step
    # but the spacial phase factor is too expensive to calculate at runtime
    sin_decay, cos_decay = get_k_sin_cos(T, h.decay.k, dx, nele)
    sin_drive, cos_drive = get_drives_sin_cos(h, dx, nele)

    phase_drive = get_drives_tracker(h)

    # Somehow helps type inference. Might be a typeinf bug or because
    # there's too many parameters
    P = typeof(p_fft!)
    PI = typeof(p_bfft!)
    SystemPropagator{H, T, N, P, PI}(h, dt, dx, nstep, nele,
                                     tmp, p_fft!, p_bfft!,
                                     E_k, (E_xg, E_xe), P_k, (P_x2g, P_x2e),
                                     sin_decay, cos_decay,
                                     sin_drive, cos_drive, phase_drive)
end

@enum AccumType AccumX AccumK

abstract AbstractAccumulator

@inline accum_init(::AbstractAccumulator, ::SystemPropagator) = nothing

function accumulate
end

# Before the iterations start the accumulator is called with
#     accum_init(accumulator, P)
# This should be used to initialize internal states (e.g. buffers).
#
# At each iteration accumulator is called with
#     accumulate(accumulator, P, i, ψ, accum_type::AccumType)
# Where
# P is the propagator
# i is the time index (1:(nstep + 1))
# ψ is the current wavefunction
# accum_type is the type of the wavefunction (X or K basis)
function propagate{H, T, N}(P::SystemPropagator{H, T, N},
                            ψ0::Matrix{Complex{T}}, # 2 x nele
                            accumulator::AbstractAccumulator)
    accum_init(accumulator, P)
    # Disable denormal values
    ccall(:jl_zero_subnormals, UInt8, (UInt8,), 1)
    eΓ4 = exp(-P.H.decay.Γ * P.dt / 4)

    ψ_norm::T = 0
    @inbounds for i in 1:P.nele
        ψ_g = ψ0[1, i]
        ψ_e = ψ0[2, i]
        P.tmp[1, i] = ψ_g
        P.tmp[2, i] = ψ_e
        ψ_norm += abs2(ψ_g) + abs2(ψ_e)
    end

    @inbounds for i in 1:(P.nstep + 1)
        ψ_scale = 1 / sqrt(ψ_norm)
        p_decay::T = 0
        for j in 1:P.nele
            ψ_g = P.tmp[1, j] * P.P_x2[1][j] * ψ_scale
            ψ_e = P.tmp[2, j] * P.P_x2[2][j] * ψ_scale
            P.tmp[1, j] = ψ_g
            P.tmp[2, j] = ψ_e
            p_decay += abs2(ψ_e)
        end
        if rand() < p_decay * P.H.decay.Γ * P.dt
            ksign = rand() > 0.5 ? 1im : -1im
            ψ_scale = 1 / sqrt(p_decay)
            for j in 1:P.nele
                ψ_e = P.tmp[2, j]
                ψ_e *= (P.cos_decay[j] + ksign * P.sin_decay[j]) * ψ_scale
                P.tmp[2, j] = 0
                P.tmp[1, j] = ψ_e
            end
            accumulate(accumulator, P, i, P.tmp, AccumX)
            P.p_fft!(P.tmp)
            # P.p_fft! * P.tmp
            ψ_scale = 1 / sqrt(T(P.nele))
            for j in 1:P.nele
                P.tmp[1, j] *= ψ_scale
                P.tmp[2, j] *= ψ_scale
            end
            accumulate(accumulator, P, i, P.tmp, AccumK)
            P.p_bfft!(P.tmp)
            # P.p_bfft! * P.tmp
            ψ_norm = T(P.nele)
            continue
        end
        accumulate(accumulator, P, i, P.tmp, AccumX)
        P.p_fft!(P.tmp)
        # P.p_fft! * P.tmp
        ψ_scale = 1 / sqrt(T(P.nele))
        for j in 1:P.nele
            p_k = P.P_k[j] * ψ_scale
            P.tmp[1, j] *= p_k
            P.tmp[2, j] *= p_k
        end
        accumulate(accumulator, P, i, P.tmp, AccumK)
        P.p_bfft!(P.tmp)
        # P.p_bfft! * P.tmp
        ψ_norm = 0
        for j in 1:P.nele
            ψ_g = P.tmp[1, j] * P.P_x2[1][j]
            ψ_e = P.tmp[2, j] * P.P_x2[2][j] * eΓ4

            # T12 = P.P_σ12[j, i - 1]
            # T11 = P.P_σ11[j, i - 1]
            # T22 = P.P_σ22[j, i - 1]
            # T21 = -conj(T12)

            # ψ_g, ψ_e = (T22 * ψ_g + T21 * ψ_e), (T11 * ψ_e + T12 * ψ_g)

            ψ_e = ψ_e * eΓ4
            ψ_norm += abs2(ψ_g) + abs2(ψ_e)
            P.tmp[2, j] = ψ_e
            P.tmp[1, j] = ψ_g
        end
    end
    ccall(:jl_zero_subnormals, UInt8, (UInt8,), 0)
end

type WaveFuncRecorder{Acc, T} <: AbstractAccumulator
    ψs::Array{Complex{T}, 3}
    WaveFuncRecorder() = new()
    WaveFuncRecorder(ψs) = new(ψs)
end

function call{H, T, Acc}(::Type{WaveFuncRecorder{Acc}},
                         P::SystemPropagator{H, T})
    return WaveFuncRecorder{Acc, T}(Array{Complex{T}}(2, P.nele, P.nstep + 1))
end

function accum_init{H, T, Acc}(r::WaveFuncRecorder{Acc, T},
                               P::SystemPropagator{H, T})
    if !isdefined(r, :ψs) || size(r.ψs) != (2, P.nele, P.nstep + 1)
        r.ψs = Array{Complex{T}}(2, P.nele, P.nstep + 1)
    end
    nothing
end

function accumulate{H, T}(r::WaveFuncRecorder{AccumX, T},
                          P::SystemPropagator{H, T}, t_i,
                          ψ::Matrix{Complex{T}}, accum_type::AccumType)
    @inbounds if accum_type == AccumX
        for j in 1:P.nele
            r.ψs[1, j, t_i] = ψ[1, j]
            r.ψs[2, j, t_i] = ψ[2, j]
        end
    end
    nothing
end

function accumulate{H, T}(r::WaveFuncRecorder{AccumK, T},
                          P::SystemPropagator{H, T}, t_i,
                          ψ::Matrix{Complex{T}}, accum_type::AccumType)
    @inbounds if accum_type == AccumK
        k_off = P.nele ÷ 2
        for j in 1:P.nele
            k = j + k_off
            k = ifelse(k > P.nele, k - P.nele, k)
            r.ψs[1, k, t_i] = ψ[1, j]
            r.ψs[2, k, t_i] = ψ[2, j]
        end
    end
    nothing
end

type EnergyRecorder{T} <: AbstractAccumulator
    Es::Vector{T}
end

function call{H, T}(::Type{EnergyRecorder}, P::SystemPropagator{H, T})
    EnergyRecorder{T}(Array{T}(P.nstep + 1))
end

function accum_init{H, T}(r::EnergyRecorder{T}, P::SystemPropagator{H, T})
    fill!(r.Es, T(0))
    nothing
end

@inline function accumulate{H, T}(r::EnergyRecorder{T},
                                  P::SystemPropagator{H, T}, t_i,
                                  ψ::Matrix{Complex{T}},
                                  accum_type::AccumType)
    @inbounds if accum_type == AccumK
        for j in 1:P.nele
            r.Es[t_i] += (abs2(ψ[1, j]) + abs2(ψ[2, j])) * P.E_k[j]
        end
    else
        for j in 1:P.nele
            # Use ground state potential for both, representing the average
            # potential energy after decay.
            r.Es[t_i] += (abs2(ψ[1, j]) + abs2(ψ[2, j])) * P.E_x[1][j]
            # r.Es[t_i] += (abs2(ψ[1, j]) * P.E_x[1][j] +
            #               abs2(ψ[2, j]) * P.E_x[2][j])
        end
    end
    nothing
end

# Time unit: μs
# Length unit: μm
# Frequency unit: MHz

# m here is actually m / ħ
m_Na = 22.98977e-3 / 6.02214129e23 / (1.0545717253362894e-34 * 1e6)
ω_g = 2π * 0.1 # f = 100kHz
ω_e = 2π * 0.1 # f = 100kHz
h_trap = HTrap(m_Na, (ω_g, ω_e))

# k, Γ
o_decay = OpticalDecay(2π * 0.589, 2π * 10.0)

# k, Ω, δ
o_drive1 = OpticalDrive(2π * 0.589, 2π * 5.0, 2π * 0.0, 100.0)

h_system = HSystem(h_trap, o_decay, (o_drive1,))

grid_size = 256
grid_space = 0.01
p_sys = SystemPropagator(h_system, 0.005, grid_space, 10000, grid_size)


function gen_ψ0(grid_size, grid_space)
    x_center = (grid_size + 1) * grid_space / 2
    ψ0 = Array{Complex128}(2, grid_size)
    sum = 0.0
    @inbounds for i in 1:grid_size
        ψ = exp(-((i * grid_space - x_center + 0.2) / 0.2)^2)
        sum += abs2(ψ)
        ψ0[1, i] = ψ
        ψ0[2, i] = 0
    end
    sum = sqrt(sum)
    @inbounds for i in 1:grid_size
        ψ0[1, i] /= sum
    end
    ψ0
end

ψ0 = gen_ψ0(grid_size, grid_space)

const wf_x_accum = WaveFuncRecorder{AccumX}(p_sys)
const wf_k_accum = WaveFuncRecorder{AccumK}(p_sys)
const wf_accum = wf_x_accum
const e_accum = EnergyRecorder(p_sys)
const _accum = e_accum

println("start")

@time propagate(p_sys, ψ0, _accum)

gc()
@time propagate(p_sys, ψ0, _accum)

# ψs = wf_accum.ψs

# img = Array{Float64}(grid_size, size(ψs, 3))

# for i in 1:size(img, 2)
#     sum = 0.0
#     @inbounds for j in 1:size(img, 1)
#         img[j, i] = abs2(ψs[1, j, i]) + abs2(ψs[2, j, i])
#         # img[j, i] = ψs[1, j, i] + ψs[2, j, i]
#         sum += img[j, i]
#     end
#     # println((i, sum))
# end

using PyPlot

# figure()
# imshow(img[:, 1:10:end])
# colorbar()

# figure()
# imshow(log((img[:, 1:10:end])))
# colorbar()

fig = figure()
plot(e_accum.Es)
ylim(0, ylim()[2] * 1.1)

show()
