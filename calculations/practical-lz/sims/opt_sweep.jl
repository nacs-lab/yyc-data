#!/usr/bin/julia -f

import TwoLevels: AbstractDrive, Drives
import TwoLevels: SeqBuilder, Builders
import TwoLevels: propagate
import TwoLevels: DummyMeasure, FullMeasure, SingleMeasure
using Optim

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
    Ω
end

function gen_seq_drive(dt, n, sinN)
    builder = SeqBuilder(dt)

    drive = SinsDrive(sinN, zero(dt), zero(dt), zero(dt))
    # Builders.start_measure(builder, FullMeasure)
    Builders.add_drive(builder, n * dt, drive)
    Builders.start_measure(builder, SingleMeasure)
    Builders.finish_measure(builder)
    seq = Builders.build(builder)
    seq, drive
end

function get_end_ext{T}(seq, drive::SinsDrive{T})
    y = T[1f0, 0]
    propagate(seq, y)
    measure = seq.measure[1].second
    abs2((measure::SingleMeasure{T}).y[2])
end

# μs, MHz
immutable OptParams{T,Seq,Dri}
    seq::Seq
    dri::Dri
    Δδ::Base.RefValue{T}
    # Maximize the average between ±δ1
    # Minimize the maximum between [δ2, δ3] and [-δ3, -δ2]
    δ1::T
    δ2::T
    δ3::T
end
function call{T}(::Type{OptParams}, dt::T, n, sinN, δ1, δ2, δ3)
    @assert 0 < δ1 < δ2 < δ3
    seq, dri = gen_seq_drive(dt, n, sinN)
    OptParams{T,typeof(seq),typeof(dri)}(seq, dri, Ref{T}(0), δ1, δ2, δ3)
end
function set_params(opt::OptParams, params)
    # params: [Δδ; cδ[sinN]; cΩ[sinN]]
    dri = opt.dri
    seq = opt.seq
    cδ = dri.cδ
    cΩ = dri.cΩ
    sinN = length(cδ)
    @assert length(cΩ) == sinN
    @assert length(params) == 2 * sinN + 1
    dri.Ω0 = 0
    dri.Ω1 = 0
    @inbounds opt.Δδ[] = params[1]
    @inbounds for i in 1:sinN
        cδ[i] = params[1 + i]
    end
    @inbounds for i in 1:sinN
        cΩ[i] = params[1 + sinN + i]
    end
end
function run_range{T}(opt::OptParams{T}, δs)
    dri = opt.dri
    seq = opt.seq
    Δδ = opt.Δδ[]
    np = length(δs)
    ext_min = one(T)
    ext_max = zero(T)
    ext_avg = zero(T)
    for δ in δs
        dri.δ0 = δ - Δδ / 2
        dri.δ1 = δ + Δδ / 2
        ext = get_end_ext(seq, dri)
        ext_min = min(ext, ext_min)
        ext_max = max(ext, ext_max)
        ext_avg += ext
    end
    ext_avg /= np
    ext_min, ext_max, ext_avg
end
function gen_scan_plot{T}(opt::OptParams{T})
    δ1 = opt.δ1
    δ2 = opt.δ2
    δ3 = opt.δ3
    dri = opt.dri
    seq = opt.seq
    Δδ = opt.Δδ[]

    np = 1000 # number of samples in total
    δs = linspace(-δ3, δ3, np)
    exts = Vector{T}(np)
    @inbounds for i in 1:np
        δ = δs[i]
        dri.δ0 = δ - Δδ / 2
        dri.δ1 = δ + Δδ / 2
        exts[i] = get_end_ext(seq, dri)
    end
    δs, exts
end

const opt_params = OptParams(2f-1, 500, 8, 0.15f-1, 1f-1, 8f-1)
