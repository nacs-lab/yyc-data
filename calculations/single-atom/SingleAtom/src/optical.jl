#!/usr/bin/julia -f

# Properties related to a optical beam

module Optical

using Compat
using ..Utils

export Drive, PhaseTracker, init_phase!, update_phase!

# Amplitude

"""
Optical drive. Parameters determines

* Phase evolution (both space and time)

* Amplitude (including polarization)
"""
immutable Drive{Amp,T} # Amp::Vec3D{Complex{T}}
    k::T
    δ::T
    ϕ0::T
    τ_θ::T
    @inline function Drive(k, δ, ϕ0, τ_θ)
        Base.typeassert(Amp, Vec3D{Complex{T}})
        new(k, δ, ϕ0, τ_θ)
    end
end

@generated get_drive_amplitude{T<:Drive}(::Type{T}) = T.parameters[1]::Vec3D

@inline get_drive_type{T}(::Vec3D{Complex{T}}) = T

@compat @generated function (::Type{Drive{Amp}}){Amp}(args...)
    @meta_expr inline
    quote
        $(Expr(:meta, :inline))
        Drive{Amp,$(get_drive_type(Amp))}(args...)
    end
end

@compat (::Type{TrigCache}){Amp,T}(drive::Drive{Amp,T}, xs) =
    TrigCache{T}(xs .* drive.k)

"""
Keep track of the temporal phase evolution of an optical drive
"""
type PhaseTracker{T}
    δ::T
    ϕ0::T
    τ_θ::T

    phase::T
    exp_t::Complex{T}
end

@compat (::Type{PhaseTracker}){Amp,T}(drive::Drive{Amp,T}) =
    PhaseTracker{T}(drive.δ, drive.ϕ0, drive.τ_θ, 0, 1)

"""
Initialize the phase tracker. This reset the phase to it's initial value,
which will be ϕ0, if it is a finite number and random otherwise. end

This needs to be done before every iteration.
"""
@inline function init_phase!{T}(track::PhaseTracker{T})
    # force inline to avoid jlcall signature
    if isfinite(track.ϕ0)
        track.phase = track.ϕ0
    else
        track.phase = rand(T) * T(2π)
    end
    track.exp_t = exp(im * track.phase)
    track
end

"""
Forward propagate the phase by dt, returns the phase and the exponential
of the phase
"""
@inline function update_phase!{T}(track::PhaseTracker{T}, dt::T)
    phase = track.phase
    if isfinite(track.τ_θ)
        δτ = dt / track.τ_θ
        phase += sqrt(δτ) * (rand(T) - T(0.5)) * π
    end
    phase = (phase - track.δ * dt) % T(2π)
    exp_t = Complex(cos(phase), sin(phase))
    track.exp_t = exp_t
    track.phase = phase
    (phase, exp_t)
end

end
