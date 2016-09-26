#!/usr/bin/julia -f

module Drives

# Abstract interface for the drive
abstract AbstractDrive{T<:AbstractFloat}

function getδ end
function getΩ end
getϕ₀{T}(::AbstractDrive{T}) = zero(T)
Base.reset(::AbstractDrive) = nothing

# Phase tracker
type PhaseTracker{T<:AbstractFloat}
    δ::T
    ϕ::T
    PhaseTracker() = new(zero(T), zero(T))
end
function Base.reset{T}(tracker::PhaseTracker{T})
    tracker.δ = zero(T)
    tracker.ϕ = zero(T)
    nothing
end

# Update the current phase and detuning
function update{T}(tracker::PhaseTracker{T}, drive::AbstractDrive{T},
                   dt::T, t::T, tlen::T, vold::T)
    t1 = muladd(T(-0.8872983346207417), dt, t)
    t2 = muladd(T(-0.5), dt, t)
    t3 = muladd(T(-0.1127016653792583), dt, t)
    δ1 = getδ(drive, t1, tlen, vold)
    δ2 = getδ(drive, t2, tlen, vold)
    δ3 = getδ(drive, t3, tlen, vold)
    δ = getδ(drive, t, tlen, vold)
    dϕ = muladd(T(0.2777777777777778), δ1 + δ3,
                T(0.4444444444444444) * δ2) * dt
    tracker.ϕ += dϕ
    tracker.δ = δ
    nothing
end

getδ{T}(tracker::PhaseTracker{T}, ::AbstractDrive{T}, ::T, ::T, ::T) = tracker.δ
getΩ{T}(::PhaseTracker{T}, drive::AbstractDrive{T}, t::T, tlen::T, vold::T) =
    getΩ(drive, t, tlen, vold)
getϕ{T}(tracker::PhaseTracker{T}) = tracker.ϕ

# Constant drive
immutable ConstDrive{T<:AbstractFloat} <: AbstractDrive{T}
    δ::T
    Ω::T
    ConstDrive(δ, Ω) = new(δ, Ω)
end
ConstDrive{T<:AbstractFloat}(δ::T, Ω::T) = ConstDrive{T}(δ, Ω)
getδ{T}(drive::ConstDrive{T}, ::T, ::T, ::T) = drive.δ
getΩ{T}(drive::ConstDrive{T}, ::T, ::T, ::T) = drive.Ω

# Some more efficient specialization for the phase tracker
@inline function update{T}(tracker::PhaseTracker{T}, drive::ConstDrive{T},
                           dt::T, ::T, ::T, ::T)
    tracker.ϕ += drive.δ * dt
    tracker.δ = drive.δ
    nothing
end
getδ{T}(::PhaseTracker{T}, drive::ConstDrive{T}, ::T, ::T, ::T) = drive.δ

immutable LinearRampDrive{T<:AbstractFloat} <: AbstractDrive{T}
    δ0::T
    δ1::T
    Ω0::T
    Ω1::T
    LinearRampDrive(δ0, δ1, Ω0, Ω1=Ω0) = new(δ0, δ1, Ω0, Ω1)
end
LinearRampDrive{T<:AbstractFloat}(δ0::T, δ1::T, Ω0::T, Ω1::T=Ω0) =
    LinearRampDrive{T}(δ0, δ1, Ω0, Ω1)
getδ{T}(drive::LinearRampDrive{T}, t::T, len::T, ::T) =
    (drive.δ0 * (len - t) + drive.δ1 * t) / len
getΩ{T}(drive::LinearRampDrive{T}, t::T, len::T, ::T) =
    (drive.Ω0 * (len - t) + drive.Ω1 * t) / len

immutable RampToDrive{T<:AbstractFloat} <: AbstractDrive{T}
    δ::T
    Ω::T
    RampToDrive(δ, Ω) = new(δ, Ω)
end
RampToDrive{T<:AbstractFloat}(δ::T, Ω::T) = RampToDrive{T}(δ, Ω)
getδ{T}(drive::RampToDrive{T}, t::T, len::T, vold::T) =
    (vold * (len - t) + drive.δ * t) / len
getΩ{T}(drive::RampToDrive{T}, t::T, len::T, vold::T) =
    (vold * (len - t) + drive.Ω * t) / len

type SinsDrive{T<:AbstractFloat} <: AbstractDrive{T}
    cδ::Vector{T}
    cΩ::Vector{T}
    δ0::T
    δ1::T
    Ω0::T
    Ω1::T
    SinsDrive(N, δ0, δ1, Ω0, Ω1=Ω0) =
        new(zeros(T, N), zeros(T, N), δ0, δ1, Ω0, Ω1)
end
SinsDrive{T<:AbstractFloat}(N, δ0::T, δ1::T, Ω0::T, Ω1::T=Ω0) =
    SinsDrive{T}(N, δ0, δ1, Ω0, Ω1)
function Drives.getδ{T}(drive::SinsDrive{T}, t::T, len::T, ::T)
    δ = (drive.δ0 * (len - t) + drive.δ1 * t) / len
    θ = t / len * π
    cδ = drive.cδ
    @inbounds for i in 1:length(drive.cδ)
        δ += cδ[i] * sin(θ * i)
    end
    δ
end
function Drives.getΩ{T}(drive::SinsDrive{T}, t::T, len::T, ::T)
    Ω = (drive.Ω0 * (len - t) + drive.Ω1 * t) / len
    θ = t / len * π
    cΩ = drive.cΩ
    @inbounds for i in 1:length(drive.cΩ)
        Ω += cΩ[i] * sin(θ * i)
    end
    abs(Ω)
end

# immutable SplineDrive{T<:AbstractFloat} <: AbstractDrive{T}
#     c::Matrix{T}
#     SplineDrive(n) = new(Matrix{T}(4, n))
# end
# function Base.setindex!(dri::SplineDrive, vals)
# end

end
