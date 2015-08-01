#!/usr/bin/julia -f

# Properties related to a optical beam

module Optical

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

@inline get_drive_type{T}(::Vec3D{Complex{T}}) = T

@generated function call{Amp}(::Type{Drive{Amp}}, args...)
    @meta_expr inline
    quote
        $(Expr(:meta, :inline))
        Drive{Amp,$(get_drive_type(Amp))}(args...)
    end
end

"""
Keep track of the temporal phase evolution of an optical drive
"""
type PhaseTracker{Amp, T}
    drive::Drive{Amp, T}

    phase::T
    exp_t::Complex{T}
end

call{Amp, T}(::Type{PhaseTracker}, drive::Drive{Amp, T}) =
    PhaseTracker{Amp, T}(drive, 0, 1)

"""
Initialize the phase tracker. This reset the phase to it's initial value,
which will be ϕ0, if it is a finite number and random otherwise. end

This needs to be done before every iteration.
"""
function init_phase!{Amp, T}(track::PhaseTracker{Amp, T})
    if isfinite(track.drive.ϕ0)
        track.phase = track.drive.ϕ0
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
function update_phase!{Amp, T}(track::PhaseTracker{Amp, T}, dt::T)
    drive = track.drive
    phase = track.phase
    if isfinite(drive.τ_θ)
        δτ = dt / drive.τ_θ
        phase += sqrt(δτ) * (rand(T) - T(0.5)) * π
    end
    phase = (phase - drive.δ * dt) % T(2π)
    exp_t = exp(im * phase)
    track.exp_t = exp_t
    track.phase = phase
    (phase, exp_t)
end

end
