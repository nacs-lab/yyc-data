#!/usr/bin/julia -f

using StructsOfArrays
using Compat

typealias ComplexSoArray{T,N} StructOfArrays{Complex{T},N,NTuple{2,Array{T,N}}}
typealias ComplexSoVector{T} ComplexSoArray{T,1}
typealias ComplexSoMatrix{T} ComplexSoArray{T,2}

@inline function Base.unsafe_copy!{T,N}(dest::ComplexSoArray{T,N},
                                        src::Array{Complex{T},N})
    len = length(src)
    BLAS.blascopy!(len, Ptr{T}(pointer(src)), 2,
                   pointer(dest.arrays[1]), 1)
    BLAS.blascopy!(len, Ptr{T}(pointer(src)) + sizeof(T), 2,
                   pointer(dest.arrays[2]), 1)
    dest
end

@inline function Base.unsafe_copy!{T,N}(dest::Array{Complex{T},N},
                                        src::ComplexSoArray{T,N})
    len = length(src)
    BLAS.blascopy!(len, pointer(src.arrays[1]), 1,
                   Ptr{T}(pointer(dest)), 2)
    BLAS.blascopy!(len, pointer(src.arrays[2]), 1,
                   Ptr{T}(pointer(dest)) + sizeof(T), 2)
    dest
end

typealias BitsT32 Union{Float32,Int32,UInt32}

@inline function Base.unsafe_copy!{T<:BitsT32,N}(dest::ComplexSoArray{T,N},
                                                 src::Array{Complex{T},N})
    len = length(src)
    dest_ptr1 = Ptr{UInt32}(pointer(dest.arrays[1]))
    dest_ptr2 = Ptr{UInt32}(pointer(dest.arrays[2]))
    src_ptr = Ptr{UInt64}(pointer(src))
    @inbounds @simd for i in 1:len
        src_v = unsafe_load(src_ptr, i)
        dest_v1 = src_v % UInt32
        dest_v2 = (src_v >> 32) % UInt32
        unsafe_store!(dest_ptr1, dest_v1, i)
        unsafe_store!(dest_ptr2, dest_v2, i)
    end
    dest
end

@inline function Base.unsafe_copy!{T<:BitsT32,N}(dest::Array{Complex{T},N},
                                                 src::ComplexSoArray{T,N})
    len = length(dest)
    src_ptr1 = Ptr{UInt32}(pointer(src.arrays[1]))
    src_ptr2 = Ptr{UInt32}(pointer(src.arrays[2]))
    dest_ptr = Ptr{UInt64}(pointer(dest))
    @inbounds @simd for i in 1:len
        src_v1 = unsafe_load(src_ptr1, i)
        src_v2 = UInt64(unsafe_load(src_ptr2, i))
        dest_v = (src_v2 << 32) | src_v1
        unsafe_store!(dest_ptr, dest_v, i)
    end
    dest
end

# 1D harmonic trap (non-magic wavelength) using Monte Carlo

##
# Hamiltonian and jump operators

immutable HTrap{T}
    m::T
    ω::NTuple{2, T}
end

@inline function get_potential(h::HTrap, x, idx)
    ω = h.ω[idx]
    sign(ω) * h.m * ω^2 * x^2 / 2
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
    exp_t::Complex{T}

    sindθ_cache::T
    cosdθ_cache::T
end

@compat (::Type{PhaseTracker}){T}(drive::OpticalDrive{T}) =
    PhaseTracker{T}(drive, T(0), T(0),
                    T(0), Complex{T}(0),
                    T(0), T(1))

function phase_tracker_init{T}(track::PhaseTracker{T})
    if isfinite(track.drive.τ_θ)
        track.phase = rand(T) * T(2π)
    else
        track.phase = 0
    end
    track.prev_t = 0
    track.total_phase = track.phase
    track.exp_t = exp(im * track.total_phase)

    track.sindθ_cache = 0
    track.cosdθ_cache = 1

    track
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
    δθ = sqrt(δt) * (rand(T) - T(0.5)) * π
    track.phase = (track.phase + δθ) % T(2π)
    track.phase
end

function phase_tracker_update{T}(track::PhaseTracker{T}, t::T, dt::T)
    phase = phase_tracker_next(track, t)
    track.total_phase = (phase - track.drive.δ * t) % T(2π)
    track.exp_t = exp(im * track.total_phase)
    θdt = track.drive.Ω * dt
    track.sindθ_cache = sin(θdt)
    track.cosdθ_cache = cos(θdt)
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

    tmp::Matrix{Complex{T}} # nele x 2
    sotmp::ComplexSoMatrix{T} # nele x 2
    p_fft!::P
    p_bfft!::PI

    # Energy
    E_k::Vector{T} # nele
    E_x::NTuple{2, Vector{T}} # nele

    # Propagator
    P_k::ComplexSoVector{T} # nele
    P_x2::NTuple{2, ComplexSoVector{T}} # nele

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

@generated function init_drives_tracker{H, T, N}(p::SystemPropagator{H, T, N})
    body = Expr(:block)
    resize!(body.args, N + 1)
    @inbounds for i in 1:N
        body.args[i] = :(phase_tracker_init(p.drive_phase[$i]))
    end
    body.args[N + 1] = nothing
    body
end

@generated function update_drives_tracker{N, T}(drive_phase::NTuple{N},
                                                t::T, dt::T)
    body = Expr(:block)
    resize!(body.args, N + 1)
    @inbounds for i in 1:N
        body.args[i] = :(phase_tracker_update(drive_phase[$i], t, dt))
    end
    body.args[N + 1] = nothing
    body
end

function SystemPropagator{H, T, N}(h::HSystem{H, T, N}, dt::T, dx::T,
                                   nstep, nele)
    tmp = Matrix{Complex{T}}(nele, 2) # nele x 2

    # FFT plan
    p_fft! = plan_fft!(tmp, 1, flags=FFTW.MEASURE)
    p_bfft! = plan_bfft!(tmp, 1, flags=FFTW.MEASURE)

    E_k = Vector{T}(nele)
    E_xg = Vector{T}(nele)
    E_xe = Vector{T}(nele)

    P_k = StructOfArrays(Complex{T}, nele)
    P_x2g = StructOfArrays(Complex{T}, nele)
    P_x2e = StructOfArrays(Complex{T}, nele)

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
                                     tmp, convert(StructOfArrays, tmp),
                                     p_fft!, p_bfft!,
                                     E_k, (E_xg, E_xe), P_k, (P_x2g, P_x2e),
                                     sin_decay, cos_decay,
                                     sin_drive, cos_drive, phase_drive)
end

@enum AccumType AccumX AccumK

abstract AbstractAccumulator

@inline accum_init(::AbstractAccumulator, ::Any) = nothing
@inline accum_finalize(::AbstractAccumulator, ::Any) = nothing

function accumulate
end

@enum DecayType DecayNone DecayLeft DecayRight DecayMiddle

# Before the iterations start the accumulator is called with
#     accum_init(accumulator, P)
# This should be used to initialize internal states (e.g. buffers).
#
# At each iteration accumulator is called with
#     accumulate(accumulator, P, i, ψ, accum_type::AccumType, decay)
# Where
# P is the propagator
# i is the time index (1:(nstep + 1))
# ψ is the current wavefunction
# accum_type is the type of the wavefunction (X or K basis)
# decay::DecayType is the type of the decay (and whether or not it happened)
function propagate{H, T, N}(P::SystemPropagator{H, T, N},
                            ψ0::Matrix{Complex{T}}, # 2 x nele
                            accumulator::AbstractAccumulator)
    # @code_llvm propagate(P, ψ0, accumulator)
    accum_init(accumulator, P)
    # Disable denormal values
    set_zero_subnormals(true)
    eΓ4 = exp(-P.H.decay.Γ * P.dt / 4)
    init_drives_tracker(P)

    dt = P.dt
    tmp = P.tmp
    sotmp = P.sotmp
    P_x2_1 = P.P_x2[1]
    P_x2_2 = P.P_x2[2]
    P_k = P.P_k
    drive_phase = P.drive_phase
    sin_drive = P.sin_drive
    cos_drive = P.cos_drive
    cos_decay = P.cos_decay
    sin_decay = P.sin_decay

    ψ_norm::T = 0
    @inbounds for i in 1:P.nele
        ψ_g = ψ0[1, i]
        ψ_e = ψ0[2, i]
        sotmp[i, 1] = ψ_g
        sotmp[i, 2] = ψ_e
        ψ_norm += abs2(ψ_g) + abs2(ψ_e)
    end

    @inbounds for i in 1:(P.nstep + 1)
        ψ_scale = 1 / sqrt(ψ_norm)
        p_decay::T = 0
        @simd for j in 1:P.nele
            ψ_g = sotmp[j, 1] * P_x2_1[j] * ψ_scale
            ψ_e = sotmp[j, 2] * P_x2_2[j] * ψ_scale
            sotmp[j, 1] = ψ_g
            sotmp[j, 2] = ψ_e
            p_decay += abs2(ψ_e)
        end
        if rand() < p_decay * P.H.decay.Γ * dt
            decay_direction = rand()
            if decay_direction < T(0.25)
                ksign = 1
                decay_type = DecayRight
            elseif decay_direction < T(0.75)
                ksign = 0
                decay_type = DecayMiddle
            else
                ksign = -1
                decay_type = DecayLeft
            end
            ψ_scale = 1 / sqrt(p_decay)
            if ksign == 0
                @simd for j in 1:P.nele
                    ψ_e = sotmp[j, 2]
                    ψ_e *= ψ_scale
                    sotmp[j, 2] = 0
                    sotmp[j, 1] = ψ_e
                end
            else
                @simd for j in 1:P.nele
                    ψ_e = sotmp[j, 2]
                    ψ_e *= Complex(cos_decay[j],
                                    ksign * sin_decay[j]) * ψ_scale
                    sotmp[j, 2] = 0
                    sotmp[j, 1] = ψ_e
                end
            end
            accumulate(accumulator, P, i, sotmp, AccumX, decay_type)
            Base.unsafe_copy!(tmp, sotmp)
            P.p_fft! * tmp
            ψ_scale = 1 / sqrt(T(P.nele))
            scale!(ψ_scale, tmp)
            accumulate(accumulator, P, i, tmp, AccumK, decay_type)
            # P.p_bfft! * tmp
            # ψ_norm = T(P.nele)
            ψ_norm = T(1)
            continue
        end
        accumulate(accumulator, P, i, sotmp, AccumX, DecayNone)
        Base.unsafe_copy!(tmp, sotmp)
        P.p_fft! * tmp
        ψ_scale = 1 / sqrt(T(P.nele))
        @simd for j in 1:P.nele
            p_k = P_k[j] * ψ_scale
            tmp[j, 1] *= p_k
            tmp[j, 2] *= p_k
        end
        accumulate(accumulator, P, i, tmp, AccumK, DecayNone)
        P.p_bfft! * tmp
        update_drives_tracker(drive_phase, (i + 1) * dt, dt)
        Base.unsafe_copy!(sotmp, tmp)
        ψ_norm = 0
        eΓ2 = (eΓ4^2)
        @simd for j in 1:P.nele
            ψ_g = sotmp[j, 1] * P_x2_1[j]
            ψ_e = sotmp[j, 2] * P_x2_2[j] * eΓ2
            sotmp[j, 1] = ψ_g
            sotmp[j, 2] = ψ_e
            ψ_norm += abs2(ψ_g) + abs2(ψ_e)
        end
        for k in 1:N
            tracker = drive_phase[k]
            sins = sin_drive[k]
            coss = cos_drive[k]
            # Hamiltonian of the spin part is
            # H_σ = Ω (cos(θ_t + θ_x) σ_x + sin(θ_t + θ_x) σ_y)

            # Propagator is
            # P_σ = exp(im H_σ Δt)
            #     = exp(im Ω (cos(θ_t + θ_x) σ_x +
            #                 sin(θ_t + θ_x) σ_y) Δt)
            #     = cos(Ω Δt) + im * (cos(θ_t + θ_x) σ_x +
            #                         sin(θ_t + θ_x) σ_y) * sin(Ω Δt)
            #     = cos(Ω Δt) + im cos(θ_t + θ_x) sin(Ω Δt) σ_x +
            #       im sin(θ_t + θ_x) sin(Ω Δt) σ_y
            #     = [cos(Ω Δt), im exp(im(θ_t + θ_x)) sin(Ω Δt)
            #        im exp(-im(θ_t + θ_x)) sin(Ω Δt), cos(Ω Δt)]

            sin_dt = tracker.sindθ_cache
            cos_dt = tracker.cosdθ_cache

            # P_σ11 = P_σ22 = cos(Ω Δt)
            # P_σ12 = im exp(im θ_t) exp(im θ_x) sin(Ω Δt)
            # P_σ21 = -P_σ12'

            exp_θ_t = tracker.exp_t
            T11 = T22 = cos_dt
            T_pre = im * exp_θ_t * sin_dt
            @simd for j in 1:P.nele
                ψ_g = sotmp[j, 1]
                ψ_e = sotmp[j, 2]
                exp_θ_x = Complex(coss[j], sins[j])

                T12 = T_pre * exp_θ_x
                T21 = -conj(T12)

                sotmp[j, 1] = T11 * ψ_g + T12 * ψ_e
                sotmp[j, 2] = T22 * ψ_e + T21 * ψ_g
            end
        end
    end
    set_zero_subnormals(false)
    accum_finalize(accumulator, P)
end

abstract MonteCarloAccumulator <: AbstractAccumulator

function propagate{H, T, N}(P::SystemPropagator{H, T, N},
                            ψ0::Matrix{Complex{T}}, # 2 x nele
                            accumulator::MonteCarloAccumulator)
    sub_accum, n = accum_init(accumulator, P)
    for i in 1:n
        propagate(P, ψ0, sub_accum)
        accumulate(accumulator, P, sub_accum)
    end
    accum_finalize(accumulator, P)
end

type WaveFuncRecorder{Acc, T} <: AbstractAccumulator
    ψs::Array{Complex{T}, 3}
    WaveFuncRecorder() = new()
    WaveFuncRecorder(ψs) = new(ψs)
end

@compat (::Type{WaveFuncRecorder{Acc}}){H,T,Acc}(P::SystemPropagator{H,T}) =
    WaveFuncRecorder{Acc, T}(Array{Complex{T}}(2, P.nele, P.nstep + 1))

function accum_init{H, T, Acc}(r::WaveFuncRecorder{Acc, T},
                               P::SystemPropagator{H, T})
    if !isdefined(r, :ψs) || size(r.ψs) != (2, P.nele, P.nstep + 1)
        r.ψs = Array{Complex{T}}(2, P.nele, P.nstep + 1)
    end
    nothing
end

function accumulate{H, T}(r::WaveFuncRecorder{AccumX, T},
                          P::SystemPropagator{H, T}, t_i,
                          ψ, accum_type::AccumType, decay)
    @inbounds if accum_type == AccumX
        for j in 1:P.nele
            r.ψs[1, j, t_i] = ψ[j, 1]
            r.ψs[2, j, t_i] = ψ[j, 2]
        end
    end
    nothing
end

function accumulate{H, T}(r::WaveFuncRecorder{AccumK, T},
                          P::SystemPropagator{H, T}, t_i,
                          ψ, accum_type::AccumType,decay)
    @inbounds if accum_type == AccumK
        k_off = P.nele ÷ 2
        for j in 1:P.nele
            k = j + k_off
            k = ifelse(k > P.nele, k - P.nele, k)
            r.ψs[1, k, t_i] = ψ[j, 1]
            r.ψs[2, k, t_i] = ψ[j, 2]
        end
    end
    nothing
end

type EnergyRecorder{T} <: AbstractAccumulator
    Es::Vector{T}
    e_thresh::T
    t_esc::T
    dt::T
    pcount::T
end

@compat (::Type{EnergyRecorder}){H,T}(P::SystemPropagator{H,T}, e_thresh) =
    EnergyRecorder{T}(Array{T}(P.nstep + 1), e_thresh, 0, P.dt, 0)

function accum_init{H, T}(r::EnergyRecorder{T}, P::SystemPropagator{H, T})
    fill!(r.Es, T(0))
    r.t_esc = 0
    r.pcount = 0
    nothing
end

@inline function accumulate{H, T}(r::EnergyRecorder{T},
                                  P::SystemPropagator{H, T}, t_i,
                                  ψ, accum_type::AccumType, decay)
    Es = r.Es
    @inbounds if accum_type == AccumK
        E_k = P.E_k
        @simd for j in 1:P.nele
            Es[t_i] += (abs2(ψ[j, 1]) + abs2(ψ[j, 2])) * E_k[j]
        end
    else
        E_x = P.E_x[1]
        @simd for j in 1:P.nele
            # Use ground state potential for both, representing the average
            # potential energy after decay.
            Es[t_i] += (abs2(ψ[j, 1]) + abs2(ψ[j, 2])) * E_x[j]
            # r.Es[t_i] += (abs2(ψ[j, 1]) * P.E_x[1][j] +
            #               abs2(ψ[j, 2]) * P.E_x[2][j])
        end
    end
    if r.t_esc == 0
        if decay != DecayNone && accum_type == AccumK
            # Don't double count...
            r.pcount += 1
        end
        if Es[t_i] > r.e_thresh
            r.t_esc = t_i * P.dt
        end
    end
    nothing
end

function accum_finalize{H, T}(r::EnergyRecorder{T}, P::SystemPropagator{H, T})
    if r.t_esc == 0
        r.t_esc = P.dt * P.nstep
    end
    nothing
end

type WaveFuncMonteCarloRecorder{Acc, T} <: MonteCarloAccumulator
    ψs2::Array{T, 3}
    sub_accum::WaveFuncRecorder{Acc, T}
    count::Int
    ncycle::Int
end

@compat (::Type{WaveFuncMonteCarloRecorder}){Acc,T}(sub_accum::WaveFuncRecorder{Acc, T}, n) =
    WaveFuncMonteCarloRecorder{Acc,T}(Array{T}(size(sub_accum.ψs)...),
                                      sub_accum, 0, n)

function accum_init{Acc, T}(r::WaveFuncMonteCarloRecorder{Acc, T}, P)
    if size(r.ψs2) != (2, P.nele, P.nstep + 1)
        r.ψs2 = zeros(T, (2, P.nele, P.nstep + 1))
    else
        fill!(r.ψs2, 0)
    end
    r.count = 0
    r.sub_accum, r.ncycle
end

function accumulate(r::WaveFuncMonteCarloRecorder,
                    P::SystemPropagator, sub_accum)
    @assert size(r.ψs2) == size(sub_accum.ψs)
    @inbounds for i in eachindex(sub_accum.ψs)
        r.ψs2[i] += abs2(sub_accum.ψs[i])
    end
    r.count += 1
    nothing
end

function accum_finalize(r::WaveFuncMonteCarloRecorder, P)
    @inbounds for i in eachindex(r.ψs2)
        r.ψs2[i] /= r.count
    end
    r.count = 1
    nothing
end

type EnergyMonteCarloRecorder{T} <: MonteCarloAccumulator
    Es::Vector{T} # Σ(E) / average
    Es2::Vector{T} # Σ(E^2) / uncertainty
    t_esc::T # Σ(t_esc) / average
    t_esc2::T # Σ(t_esc^2) / uncertainty
    pcount::T # Σ(pcount) / average
    pcount2::T # Σ(pcount^2) / uncertainty
    sub_accum::EnergyRecorder{T}
    count::Int
    ncycle::Int
end

@compat (::Type{EnergyMonteCarloRecorder}){T}(sub_accum::EnergyRecorder{T}, n) =
    EnergyMonteCarloRecorder{T}(Array{T}(size(sub_accum.Es)),
                                Array{T}(size(sub_accum.Es)),
                                0, 0, 0, 0, sub_accum, -1, n)

function accum_init{T}(r::EnergyMonteCarloRecorder{T}, P)
    if length(r.Es2) != P.nstep + 1
        r.Es = zeros(T, P.nstep + 1)
        r.Es2 = zeros(T, P.nstep + 1)
    else
        fill!(r.Es, 0)
        fill!(r.Es2, 0)
    end
    r.t_esc = 0
    r.t_esc2 = 0
    r.pcount = 0
    r.pcount2 = 0
    r.count = 0
    r.sub_accum, r.ncycle
end

function accumulate(r::EnergyMonteCarloRecorder,
                    P::SystemPropagator, sub_accum)
    @assert r.count >= 0
    @assert size(r.Es) == size(sub_accum.Es)
    sub_Es = sub_accum.Es
    r_Es = r.Es
    r_Es2 = r.Es2
    @inbounds @simd for i in eachindex(sub_Es)
        Es = sub_Es[i]
        r_Es[i] += Es
        r_Es2[i] += Es^2
    end
    r.t_esc += sub_accum.t_esc
    r.t_esc2 += sub_accum.t_esc^2
    r.pcount += sub_accum.pcount
    r.pcount2 += sub_accum.pcount^2
    r.count += 1
    nothing
end

function sum2average(s, s2, count)
    avg = s / count
    avg2 = s2 / count
    std = (avg2 - avg^2) / (count - 1)
    # rounding errors can make small std smaller than zero
    unc = std <= 0 ? zero(std) : sqrt(std)
    avg, unc
end

function accum_finalize(r::EnergyMonteCarloRecorder, P)
    @assert r.count >= 0
    @assert size(r.Es) == size(r.Es2)
    @inbounds for i in eachindex(r.Es)
        r.Es[i], r.Es2[i] = sum2average(r.Es[i], r.Es2[i], r.count)
    end
    r.t_esc, r.t_esc2 = sum2average(r.t_esc, r.t_esc2, r.count)
    r.pcount, r.pcount2 = sum2average(r.pcount, r.pcount2, r.count)
    r.count = -1
    nothing
end

@compat (::Type{MonteCarloAccumulator})(sub_accum::EnergyRecorder, n) =
    EnergyMonteCarloRecorder(sub_accum, n)

@compat (::Type{MonteCarloAccumulator})(sub_accum::WaveFuncRecorder, n) =
    WaveFuncMonteCarloRecorder(sub_accum, n)

function plot_accum_img(img::Matrix{Float64})
    xsize, ysize = size(img)

    if xsize > ysize * 3
        xscale = xsize ÷ (ysize * 2)
        img = img[1:xscale:end, :]
    elseif ysize > xsize * 3
        yscale = ysize ÷ (xsize * 2)
        img = img[:, 1:yscale:end]
    end

    # figure()
    imshow(img)
    colorbar()

    # figure()
    # imshow(log(img))
    # colorbar()

    nothing
end

function plot_accum(accum::WaveFuncRecorder)
    ψs = accum.ψs

    img = Array{Float64}(size(ψs, 2, 3))

    for i in 1:size(img, 2)
        @inbounds for j in 1:size(img, 1)
            img[j, i] = abs2(ψs[1, j, i]) + abs2(ψs[2, j, i])
        end
    end
    plot_accum_img(img)
end

function plot_accum(accum::WaveFuncMonteCarloRecorder)
    ψs = accum.ψs2

    img = Array{Float64}(size(ψs, 2, 3))

    for i in 1:size(img, 2)
        @inbounds for j in 1:size(img, 1)
            img[j, i] = ψs[1, j, i] + ψs[2, j, i]
        end
    end
    plot_accum_img(img)
end

function plot_accum(accum::EnergyRecorder)
    # figure()
    plot(accum.Es)
    grid()
    ylim(0, ylim()[2] * 1.1)
end

function plot_accum(accum::EnergyMonteCarloRecorder)
    # figure()
    @printf("Escape time: %.2f±%.2f\n", accum.t_esc, accum.t_esc2)
    @printf("Photon Emitted: %.2f±%.2f\n", accum.pcount, accum.pcount2)
    Es = accum.Es
    Es2 = accum.Es2
    len = length(Es)
    step_len = max(len ÷ 10000, 1)
    plot_len = len ÷ step_len
    ts = (1:plot_len) * (accum.sub_accum.dt * step_len)
    avg_Es = Vector{eltype(Es)}(plot_len)
    avg_Es2 = Vector{eltype(Es2)}(plot_len)
    @inbounds for i in 1:plot_len
        es = zero(eltype(Es))
        es2 = zero(eltype(Es2))
        @simd for j in 1:step_len
            es += Es[(i - 1) * step_len + j]
            es2 += abs2(Es2[(i - 1) * step_len + j])
        end
        avg_Es[i] = es / step_len
        avg_Es2[i] = sqrt(es2 / step_len)
    end
    errorbar(ts, avg_Es, avg_Es2)
    axvline(x=accum.t_esc, color="r")
    axvline(x=accum.t_esc - accum.t_esc2, color="b", linestyle="--")
    axvline(x=accum.t_esc + accum.t_esc2, color="b", linestyle="--")
    grid()
end

nothing
