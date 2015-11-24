#!/usr/bin/julia -f

module Sequences

import ..TwoLevels: AbstractDrive, getδ, getΩ, getϕ₀
import ..TwoLevels: PhaseTracker, update, getϕ
import ..TwoLevels: AbstractMeasure, Measures

type DriveTracker{T}
    δ::T
    Ω::T
end
function Base.reset{T}(tracker::DriveTracker{T})
    tracker.δ = zero(T)
    tracker.Ω = zero(T)
    nothing
end

immutable Sequence{T,Ds<:Tuple,M<:AbstractMeasure}
    y::Vector{T}
    dt::T
    nsteps::Int
    drives::Ds
    measure::M
    phase_tracker::PhaseTracker{T}
    drive_tracker::DriveTracker{T}
    function Sequence{N}(dt, nsteps,
                         drives::NTuple{N,Pair{Int}}, measure::M)
        Base.typeassert(drives, Ds)
        y = Vector{T}(2)
        phase_tracker = PhaseTracker{T}()
        drive_tracker = DriveTracker{T}(zero(T), zero(T))
        new(y, dt, nsteps, drives, measure, phase_tracker, drive_tracker)
    end
end

function propagate!(seq::Sequence, y0, ϕ₀=0f0)
    reset(seq)
    seq.y[:] = y0
    _propagate!(seq.y, seq.drives, seq.dt, seq.measure,
                seq.phase_tracker, seq.drive_tracker, ϕ₀)
end

function Base.reset(seq::Sequence)
    fill!(seq.y, 0)
    reset_drives!(seq.drives)
    reset(seq.measure)
    reset(seq.phase_tracker)
    reset(seq.drive_tracker)
    nothing
end

@generated function reset_drives{N}(drives::NTuple{N,Pair{Int}})
    body = quote
    end
    for i in 1:N
        push!(body.args, :(reset(drives[$i].second)))
    end
    push!(body.args, nothing)
    body
end

@generated function _propagate!{N}(y, drives::NTuple{N,Pair{Int}},
                                   dt, measure::AbstractMeasure,
                                   phase_tracker, drive_tracker, ϕ₀)
    body = quote
        idx_offset = 1
        Measure.snapshot(measure, y, 1)
        phase_tracker.ϕ = ϕ₀
    end
    for i in 1:N
        ex = quote
            let
                nsteps = drives[$i].first
                drive = drives[$i].second
                propagate_step(y, drive, dt, nsteps, measure, phase_tracker,
                               drive_tracker, idx_offset)
                idx_offset += nsteps
            end
        end
        push!(body.args, ex)
    end
    push!(body.args, nothing)
    body
end

function propagate_step{T}(y::Vector{T}, drive, dt, nstep, measure,
                           phase_tracker, drive_tracker, idx_offset)
    oldδ = drive_tracker.δ
    oldΩ = drive_tracker.Ω
    phase_tracker.ϕ += getϕ₀(drive)
    tlen = nsteps * dt
    t_offset = idx_offset * dt
    @inbounds for i in 1:nsteps
        t = dt * (i - 0.5)
        update(phase_tracker, drive, ifelse(i == 1, dt / 2, dt), t, tlen, oldδ)
        # original parameters
        δ = getδ(phase_tracker, drive, t, tlen, oldδ)::T
        Ω = getΩ(phase_tracker, drive, t, tlen, oldΩ)::T
        ϕ = getϕ(phase_tracker)::T

        # derived values
        δ′ = δ
        Ω′ = √(Ω^2 + δ′^2 / 4)
        δt_2 = δ′ * dt / 2
        Ωt = Ω′ * dt

        @fastmath cosΩ = cos(Ωt)
        @fastmath sinΩ = sin(Ωt)
        @fastmath expδ = Complex(cos(δt_2), sin(δt_2))
        @fastmath expϕ = Complex(cos(ϕ), sin(ϕ))

        y_1 = y[1]
        y_2 = y[2]

        T22 = expδ * Complex(cosΩ, -δ′ / (2Ω′) * sinΩ)
        T11 = conj(T22)
        T21′ = Ω / Ω′ * expϕ * expδ * sinΩ
        T21 = Complex(imag(T21′), -real(T21′))
        T12 = Complex(-imag(T21′), -real(T21′))

        y_1′ = muladd(T12, y_2, T11 * y_1)
        y_2′ = muladd(T21, y_1, T22 * y_2)
        y_scale = 1 / sqrt(abs2(y_1) + abs2(y_2))
        y[1] = y_1′ * y_scale
        y[2] = y_2′ * y_scale
        Measures.snapshot(measure, y, i + idx_offset, t_offset + t + dt / 2)
    end
    update(phase_tracker, drive, dt / 2, tlen, tlen, oldδ)
    drive_tracker.δ = getδ(phase_tracker, drive, tlen, tlen, oldδ)::T
    drive_tracker.Ω = getΩ(phase_tracker, drive, tlen, tlen, oldΩ)::T
    nothing
end

end
