#!/usr/bin/julia -f

# Abstract interface
abstract AbstractDrive

function update_dt end
function get_detuning end
function get_rabi end
function get_phase end

abstract AbstractMeasure

function measure_snapshot end

abstract AbstractPhaseTracker

type PhaseTracker <: AbstractPhaseTracker
    ϕ::Float32
end

@inline get_phase(tracker::PhaseTracker) = tracker.ϕ
@inline update_phase(tracker::PhaseTracker, ϕ) = (tracker.ϕ = ϕ; nothing)

immutable DummyPhaseTracker <: AbstractPhaseTracker
    ϕ::Float32
end

@inline get_phase(tracker::DummyPhaseTracker) = tracker.ϕ
@inline update_phase(::DummyPhaseTracker, ϕ) = nothing

function internal_propagate!(y0::Vector, drive::AbstractDrive, dt, nsteps,
                             measure::AbstractMeasure,
                             tracker::AbstractPhaseTracker)
    ϕ = ϕ₀ = get_phase(tracker)
    @inbounds for i in 1:nsteps
        update_dt(drive, i == 1 ? dt / 2 : dt)
        # original parameters
        δ = get_detuning(drive)::Real
        Ω = get_rabi(drive)::Real
        ϕ = get_phase(drive)::Real + ϕ₀

        # derived values
        δ′ = δ
        Ω′ = √(Ω^2 + δ′^2 / 4)
        δt_2 = δ′ * dt / 2
        Ωt = Ω′ * dt

        @fastmath cosΩ = cos(Ωt)
        @fastmath sinΩ = sin(Ωt)
        @fastmath expδ = Complex(cos(δt_2), sin(δt_2))
        @fastmath expϕ = Complex(cos(ϕ), sin(ϕ))

        y_1 = y0[1]
        y_2 = y0[2]

        T22 = expδ * Complex(cosΩ, -δ′ / (2Ω′) * sinΩ)
        T11 = conj(T22)
        T21′ = Ω / Ω′ * expϕ * expδ * sinΩ
        T21 = Complex(imag(T21′), -real(T21′))
        T12 = Complex(-imag(T21′), -real(T21′))

        y_1′ = muladd(T12, y_2, T11 * y_1)
        y_2′ = muladd(T21, y_1, T22 * y_2)
        y_scale = 1 / sqrt(abs2(y_1) + abs2(y_2))
        y0[1] = y_1′ * y_scale
        y0[2] = y_2′ * y_scale
        measure_snapshot(measure, y0, i + 1)
    end
    update_dt(drive, dt / 2)
    update_phase(tracker, ϕ)
    y0
end

# Propagate function
function propagate!(y0::Vector, drive::AbstractDrive, dt, nsteps,
                    measure::AbstractMeasure, ϕ₀=0f0)
    @assert size(y0) == (2,)
    measure_snapshot(measure, y0, 1)
    internal_propagate!(y0, drive, dt, nsteps, measure, DummyPhaseTracker(ϕ₀))
end

immutable OffsetMeasure{T<:AbstractMeasure} <: AbstractMeasure
    measure::T
    offset::typeof(Ref(0))
end

@inline function measure_snapshot(measure::OffsetMeasure, y, idx)
    measure_snapshot(measure.measure, y, idx + measure.offset[])
end

@generated function propagate!{N}(y0::Vector, drives::NTuple{N,Pair{Int}},
                                  dt, measure::AbstractMeasure, ϕ₀=0f0)
    body = quote
        @assert size(y0) == (2,)
        offset_measure = OffsetMeasure(measure, Ref(0))
        tracker = PhaseTracker(ϕ₀)
        measure_snapshot(measure, y0, 1)
    end
    for i in 1:N
        ex = quote
            let
                nsteps = drives[$i].first
                drive = drives[$i].second
                internal_propagate!(y0, drive, dt, nsteps,
                                    offset_measure, tracker)
                offset_measure.offset[] += nsteps
            end
        end
        push!(body.args, ex)
    end
    push!(body.args, :y0)
    body
end

# Constant drive
immutable ConstDrive <: AbstractDrive
    δ::Float32
    Ω::Float32
    ϕ::typeof(Ref(1f0))
    ConstDrive(δ, Ω) = new(δ, Ω, Ref(0f0))
end

@inline update_dt(drive::ConstDrive, dt) = drive.ϕ[] += dt * drive.δ
@inline get_detuning(drive::ConstDrive) = drive.δ
@inline get_rabi(drive::ConstDrive) = drive.Ω
@inline get_phase(drive::ConstDrive) = drive.ϕ[]

# LZ drive
immutable LZDrive <: AbstractDrive
    δ0::Float32
    dδ::Float32
    Ω::Float32
    t::typeof(Ref(1f0))
    LZDrive(δ0, dδ, Ω) = new(δ0, dδ, Ω, Ref(0f0))
end

@inline update_dt(drive::LZDrive, dt) = drive.t[] += dt
@inline get_detuning(drive::LZDrive) = muladd(drive.dδ, drive.t[], drive.δ0)
@inline get_rabi(drive::LZDrive) = drive.Ω
@inline function get_phase(drive::LZDrive)
    t = drive.t[]
    t * muladd(drive.dδ / 2, t, drive.δ0)
end

# Full measure

immutable FullMeasure <: AbstractMeasure
    ys::Matrix{Float32}
    FullMeasure(nsteps) = new(Matrix{Float32}(2, nsteps + 1))
end

@inline function measure_snapshot(measure::FullMeasure, y, idx)
    @inbounds measure.ys[1, idx] = abs2(y[1])
    @inbounds measure.ys[2, idx] = abs2(y[2])
    nothing
end

immutable DummyMeasure <: AbstractMeasure
end

@inline function measure_snapshot(::DummyMeasure, y, idx)
    nothing
end

function main()
    y0 = Complex64[1, 0]
    lz_nsteps = 5000
    measure = FullMeasure(lz_nsteps * 2)
    drive = LZDrive(-10f0, 2f-2, 0.3f-1)
    @time propagate!(y0, drive, 1f-1, lz_nsteps * 2, measure)
    measure.ys
end

function final_value(f_c, δf, t, Ω, nsteps)
    f_start = f_c - δf / 2
    df = δf / t
    drive = LZDrive(f_start, df, Ω)
    y0 = Complex64[1, 0]
    measure = DummyMeasure()
    dt = t / nsteps
    propagate!(y0, drive, dt, nsteps, measure)
    abs2(y0[2])
end

function gen_curve()
    Ω = 3f0
    δf = 10f0
    t = 1f3
    nsteps = 10_000

    f_cs = linspace(-40f0, 0f0, 4000)
    f_cs, [final_value(f_c, δf, t, Ω, nsteps) for f_c in f_cs]
end

@time const f_cs, finals = gen_curve()

using PyPlot

figure()
plot(f_cs, finals)
xlabel("Final frequency")
ylabel("Excited state populateion")
grid()
# savefig("normal_lz.png")
show()