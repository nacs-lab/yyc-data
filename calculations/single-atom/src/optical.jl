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
        typeassert(Amp, Amplitude3D{T})
        new(k, δ, ϕ0, τ_θ)
    end
end

"""
Keep track of the phase (and other parameters) that describes the evolution
and time dependent interaction between an optical drive an a transition
"""
type PhaseTracker{Amp, T}
    drive::Drive{Amp, T}
    Ω::T # Product of the amplitude and the dipole moment

    phase::T
    exp_t::Complex{T}

    sindθ::T
    cosdθ::T
end

call{Amp, T}(::Type{PhaseTracker}, drive::Drive{Amp, T}, Ω) =
    PhaseTracker{Amp, T}(drive, Ω, 0, 0, 0, 1)

function init_phase{Amp, T}(track::PhaseTracker{Amp, T})
    if isfinite(track.drive.ϕ0)
        track.phase = track.drive.ϕ0
    else
        track.phase = rand(T) * T(2π)
    end
    track.exp_t = exp(im * track.phase)

    track.sindθ = 0
    track.cosdθ = 1

    track
end

function update_phase{Amp, T}(track::PhaseTracker{Amp, T}, dt::T)
    if isfinite(track.drive.τ_θ)
        δτ = dt / track.drive.τ_θ
        track.phase += sqrt(δτ) * (rand(T) - T(0.5)) * π
    end
    track.phase = (track.phase - track.drive.δ * dt) % T(2π)
    track.exp_t = exp(im * track.total)

    θdt = tracker.Ω * dt
    tracker.sindθ = sin(θdt)
    tracker.cosdθ = cos(θdt)
end

end
