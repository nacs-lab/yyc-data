#!/usr/bin/julia -f

module Optical

###
# Optical drive

@enum Polarization Pol_σplus Pol_σminus Pol_π

"""
Type of an optical dipole transition
"""
Polarization

"""
Amplitude of an optical wave
"""
abstract Amplitude

Base.abs(amp::Amplitude) = sqrt(abs2(amp))

"""
Amplitude vector in 3D
"""
immutable Amplitude3D{T<:Real} <: Amplitude
    x::Complex{T}
    y::Complex{T}
    z::Complex{T}
    @inline Amplitude3D(x, y, z) = new(x, y, z)
    @inline function Amplitude3D(pol::Polarization,
                                 amp::Union{T, Complex{T}}=T(1))
        if pol == Pol_π
            return new(0, 0, amp)
        end
        amp2 = sqrt(T(0.5)) * amp
        if pol == Pol_σplus
            new(amp2, im * amp2, 0)
        else
            new(amp2, -im * amp2, 0)
        end
    end
end

@inline call{T<:Real}(::Type{Amplitude3D}, pol::Polarization, amp::T) =
    Amplitude3D{T}(pol, amp)

@inline call{T<:Real}(::Type{Amplitude3D}, pol::Polarization,
                      amp::Complex{T}) = Amplitude3D{T}(pol, amp)

Base.abs2(amp::Amplitude3D) = abs2(amp.x) + abs2(amp.y) + abs2(amp.z)

immutable LinearPol
end
immutable CircularPol{Dir}
end

const linearPol = LinearPol()
const rightCirPol = CircularPol{1}()
const leftCirPol = CircularPol{-1}()

immutable Amplitude2D{T<:Real} <: Amplitude
    x::Complex{T}
    y::Complex{T}
    @inline Amplitude2D(x, y) = new(x, y)
    @inline Amplitude2D(::LinearPol, amp, angl=T(0)) =
        new(amp * cos(angl), amp * sin(angl))
    @inline function Amplitude2D{Dir}(::CircularPol{Dir}, amp)
        amp2 = sqrt(T(0.5)) * amp
        new(amp2, Dir * im * amp2)
    end
end

@inline call{T}(::Type{Amplitude2D}, p::LinearPol, amp::T, angl=T(0)) =
    Amplitude2D{T}(p, amp, angl)

@inline call{T}(::Type{Amplitude2D}, p::LinearPol, amp::Complex{T},
                angl=T(0)) = Amplitude2D{T}(p, amp, angl)

@inline call{T}(::Type{Amplitude2D}, p::CircularPol, amp::T) =
    Amplitude2D{T}(p, amp)

@inline call{T}(::Type{Amplitude2D}, p::CircularPol, amp::Complex{T}) =
    Amplitude2D{T}(p, amp)

Base.abs2(amp::Amplitude2D) = abs2(amp.x) + abs2(amp.y)

"""
Optical drive. Parameters determines

* Phase evolution (both space and time)

* Amplitude (including polarization)

"""
immutable Drive{Amp, T} # Amp::Amplitude3D
    k::T
    δ::T
    ϕ0::T
    τ_θ::T
    function Drive(k, δ, ϕ0, τ_θ)
        new(k, δ, ϕ0, τ_θ)
    end
end

type PhaseTracker{T}
    drive::OpticalDrive{T}
    phase::T
    prev_t::T

    total_phase::T
    exp_t::Complex{T}

    dθ_cached::T
    sindθ_cache::T
    cosdθ_cache::T
end

function call{T}(::Type{PhaseTracker}, drive::OpticalDrive{T})
    PhaseTracker{T}(drive, T(0), T(0),
                    T(0), Complex{T}(0),
                    T(0), T(0), T(1))
end

function phase_tracker_init{T}(track::PhaseTracker{T})
    if isfinite(track.drive.τ_θ)
        track.phase = rand(T) * T(2π)
    else
        track.phase = 0
    end
    track.prev_t = 0
    track.total_phase = track.phase
    track.exp_t = exp(im * track.total_phase)

    track.dθ_cached = 0
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

function phase_tracker_update{T}(track::PhaseTracker{T}, t::T)
    phase = phase_tracker_next(track, t)
    track.total_phase = (phase - track.drive.δ * t) % T(2π)
    track.exp_t = exp(im * track.total_phase)
    nothing
end

end
