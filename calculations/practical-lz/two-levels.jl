#!/usr/bin/julia -f

# Abstract interface
abstract AbstractDrive

function update_dt end
function get_detuning end
function get_rabi end
function get_phase end

abstract AbstractMeasure

function measure_snapshot end

# Propagate function
function propagate!(y0::Vector, drive::AbstractDrive, dt, nsteps,
                    measure::AbstractMeasure)
    # @assert size(y0) == (2,)
    measure_snapshot(measure, y0, 1)
    @inbounds for i in 1:nsteps
        # Update drive state
        update_dt(drive, dt)

        # original parameters
        δ = get_detuning(drive)::Real
        Ω = get_rabi(drive)::Real
        ϕ₀ = get_phase(drive)::Real

        # derived values
        δ′ = δ
        Ω′ = √(Ω^2 + δ′^2 / 4)
        δt_2 = δ′ * dt / 2
        Ωt = Ω′ * dt

        @fastmath cosΩ = cos(Ωt)
        @fastmath sinΩ = sin(Ωt)
        @fastmath expδ = Complex(cos(δt_2), sin(δt_2))
        @fastmath expϕ = Complex(cos(ϕ₀), sin(ϕ₀))

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
    y0
end

# Constant drive
immutable ConstDrive <: AbstractDrive
    δ::Float32
    Ω::Float32
    ϕ::typeof(Ref(1f0))
    ConstDrive(δ, Ω, ϕ₀=0f0) = new(δ, Ω, Ref(Float32(ϕ₀)))
end

@inline update_dt(drive::ConstDrive, dt) = drive.ϕ[] += dt * drive.δ
@inline get_detuning(drive::ConstDrive) = drive.δ
@inline get_rabi(drive::ConstDrive) = drive.Ω
@inline get_phase(drive::ConstDrive) = drive.ϕ[]

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

function main()
    nsteps = 1000
    y0 = Complex64[1, 0]
    drive = ConstDrive(5f0, 0.5f0)
    measure = FullMeasure(nsteps)
    @code_llvm propagate!(y0, drive, 1f-2, nsteps, measure)
    @time propagate!(y0, drive, 1f-2, nsteps, measure)
    measure.ys
end

const ys = main()
println(size(Float32[ys[2, i] for i in 1:size(ys, 2)]))

using PyPlot

figure()
plot(Float32[ys[1, i] for i in 1:size(ys, 2)], "g-")
plot(Float32[ys[2, i] for i in 1:size(ys, 2)], "r-")
show()
