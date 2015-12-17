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

function get_end_ext{T}(seq, drive::AbstractDrive{T})
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
function gen_drive_plot{T}(opt::OptParams{T})
    dri = opt.dri
    seq = opt.seq
    dt = seq.dt
    nsteps = seq.nsteps
    Δδ = opt.Δδ[]
    dri.δ0 = -Δδ / 2
    dri.δ1 = +Δδ / 2
    tlen = nsteps * dt
    δs = Vector{T}(nsteps + 1)
    Ωs = Vector{T}(nsteps + 1)

    @inbounds for i in 0:nsteps
        t = i * dt
        δs[i + 1] = Drives.getδ(dri, t, tlen, dri.δ0)
        Ωs[i + 1] = Drives.getΩ(dri, t, tlen, T(0))
    end
    (0:nsteps) * dt, δs, Ωs
end

const opt_params = OptParams(1f-1, 1000, 8, 0.2f-1, 2f-1, 8f-1)

function run_scan{T}(opt::OptParams{T})
    δ1 = opt.δ1
    δ2 = opt.δ2
    δ3 = opt.δ3
    np = 20 # number of samples in each interval
    mid_min, mid_max, mid_avg = run_range(opt, linspace(-δ1, δ1, np))
    side1_min, side1_max, side1_avg = run_range(opt, linspace(-δ3, -δ2, np))
    side2_min, side2_max, side2_avg = run_range(opt, linspace(δ2, δ3, np))
    side_min = min(side1_min, side2_min)
    side_max = max(side1_max, side2_max)
    side_avg = (side1_avg + side2_avg) / 2
    # side_min = side1_min
    # side_max = side1_max
    # side_avg = side1_avg
    (1 - mid_min) * 2 + side_max * 10 + side_avg / mid_avg * 100
end
function optim_model(params)
    # println(params)
    set_params(opt_params, params)
    res = run_scan(opt_params)
    println(res)
    Float64(res)
end

# params = [-0.11502749674683736,-0.06119036359289394,-0.034605598077178,-0.0032243129592623805,-0.0032738954954711,-0.045679977066940186,-0.0053466989997855464,-0.02194240530902273,-0.019574094058389495,0.07771069849593015,0.05225714565347806,-0.03454123505944487,-0.019839962213737725,0.00942251576695764,0.01464587601520262,-0.0013286541152173184,-0.003950156109605001]
# params = [-0.13519760613828735,-0.029440296545190148,-0.04631800963527856,-0.0002234795251961228,-0.022804659213110864,-0.03944347612559795,0.015634510300549048,-0.030649756701684138,0.008807375702501353,0.0850712530556744,0.05319441082213134,-0.03985701020388253,-0.023133794733388856,0.006171758884888809,0.008483757092033868,-0.0061279458193116846,-0.0013327518864845545]
# params = [-0.13013139682353347,0.006312375130221837,-0.031956048682332046,0.0006946369470619224,-0.0266169325552198,0.0016928797961989386,0.002108009418659553,-0.001415131715743528,-0.0022987348278135575,0.08751898962694254,0.043420080793892514,-0.04797923142566892,-0.025577461711200008,0.0101402353410935,0.006278498196434462,-0.0010224436342372658,-0.0012058875103421194]
params = [-0.13165780584535017,0.002458924776874482,-0.03886563746693689,-0.0004594378756653323,-0.016217073839701675,0.0018728848350684562,-0.006420467959789269,0.000181854329494213,-0.0030732111122789467,0.08762561850840865,0.03967701306201832,-0.04854675615963805,-0.023358415966911707,0.009540338266152673,0.007151635718551654,4.0941846022623444e-5,-0.0022496656919625334]
do_opt = true
# do_opt = false

if do_opt
    # opt_res = optimize(optim_model, params, method=:cg,
    #                    grtol=1e-6, ftol=1e-15, xtol=1e-15)
    opt_res = optimize(optim_model, params, method=:bfgs,
                       grtol=1e-6, ftol=1e-15, xtol=1e-15)
    # opt_res = optimize(optim_model, params, method=:l_bfgs,
    #                    grtol=1e-4, ftol=1e-8, xtol=1e-8)

    @show opt_res
    params = opt_res.minimum
    println(params)
end

# params: [Ω0; Ω1; Δδ; cδ[sinN]; cΩ[sinN]]
set_params(opt_params, params)
println(run_scan(opt_params))
x, y = gen_scan_plot(opt_params)
ts, δs, Ωs = gen_drive_plot(opt_params)
using PyPlot

figure()
plot(x, y, "r")
axvline(opt_params.δ1, color="b")
axvline(-opt_params.δ1, color="b")
axvline(opt_params.δ2, color="g")
axvline(-opt_params.δ2, color="g")
axvline(opt_params.δ3, color="y")
axvline(-opt_params.δ3, color="y")
grid()

figure()
plot(ts, δs, "b", label="\$\\delta\$")
plot(ts, Ωs, "g", label="\$\\Omega\$")
grid()

show()
