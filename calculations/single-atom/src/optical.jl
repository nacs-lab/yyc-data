#!/usr/bin/julia -f

# Properties related to a optical beam

module Optical

using ..Utils

export Amplitude3D, Drive, PhaseTracker, init_phase, update_phase

# Amplitude

"""
3D amplitude of an optical beam
"""
immutable Amplitude3D{T<:Real}
    x::Complex{T}
    y::Complex{T}
    z::Complex{T}
end

@inline Base.abs2(amp::Amplitude3D) = abs2(amp.x) + abs2(amp.y) + abs2(amp.z)
@inline Base.abs(amp::Amplitude3D) = sqrt(abs2(amp))

"""
Optical drive. Parameters determines

* Phase evolution (both space and time)

* Amplitude (including polarization)
"""
immutable Drive{Amp,T} # Amp::Amplitude3D
    k::T
    δ::T
    ϕ0::T
    τ_θ::T
    @inline function Drive(k, δ, ϕ0, τ_θ)
        typeassert(Amp, Amplitude3D{T})
        new(k, δ, ϕ0, τ_θ)
    end
end

@generated function call{Amp}(::Type{Drive{Amp}}, k, δ, ϕ0, τ_θ)
    @meta_expr(inline)
    amp_type = typeof(Amp::Amplitude3D)
    T = amp_type.parameters[1]
    quote
        $(Expr(:meta, :inline))
        Drive{Amp,$T}(k, δ, ψ0, τ_θ)
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
    PhaseTracker{Amp, T}(drive, 0, 0)

function init_phase{Amp, T}(track::PhaseTracker{Amp, T})
    if isfinite(track.drive.ϕ0)
        track.phase = track.drive.ϕ0
    else
        track.phase = rand(T) * T(2π)
    end
    track.exp_t = exp(im * track.phase)
    track
end

function update_phase{Amp, T}(track::PhaseTracker{Amp, T}, dt::T)
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
