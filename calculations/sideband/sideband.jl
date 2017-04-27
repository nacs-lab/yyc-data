#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 15,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

import NaCsSim: Setup, System
import NaCsCalc: Trap
using NaCsCalc.Utils: interactive
using NaCsCalc.Atomic: all_scatter_D

function maybe_save(name)
    if !interactive()
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

function plot_ground_state(params, res)
    figure()
    gp = [r.a for r in res]
    gp_unc = [r.s for r in res]
    errorbar(params, gp, gp_unc)
    title("Ground state probability")
    plot_hook()
    grid()
end

function plot_total(params, res)
    figure()
    total = [1 - t.a for t in res]
    total_unc = [t.s for t in res]
    errorbar(params, total, total_unc)
    title("Total loss")
    ylim([0, ylim()[2]])
    plot_hook()
    grid()
end

function plot_nbars(params, res)
    figure()
    nbarx_res = (nr[1] for nr in res)
    nbary_res = (nr[2] for nr in res)
    nbarz_res = (nr[3] for nr in res)

    nbarx = [n.a for n in nbarx_res]
    nbarx_unc = [n.s for n in nbarx_res]
    nbary = [n.a for n in nbary_res]
    nbary_unc = [n.s for n in nbary_res]
    nbarz = [n.a for n in nbarz_res]
    nbarz_unc = [n.s for n in nbarz_res]
    if max(maximum(nbarx), maximum(nbary), maximum(nbarz)) > 1
        gca()[:set_yscale]("log")
    end
    errorbar(params, nbarx, nbarx_unc, label="X")
    errorbar(params, nbary, nbary_unc, label="Y")
    errorbar(params, nbarz, nbarz_unc, label="Z")
    legend()
    title("\$\\bar n\$")
    plot_hook()
    grid()
end

plot_hf(params, res) = plot_hf(params, collect(res))

function plot_hf{T<:Tuple}(params, res::Vector{T})
    figure()
    for i in 1:nfields(T)
        hf = [r[i].a for r in res]
        hf_unc = [r[i].s for r in res]
        errorbar(params, hf, hf_unc, label="$i")
    end
    legend()
    ylim([0, 1])
    title("Hyperfine")
    plot_hook()
    grid()
end

function plot_result(params, res, prefix)
    plot_ground_state(params, (r[2] for r in res))
    maybe_save(joinpath(prefix, "ground.svg"))

    plot_total(params, (r[1][2] for r in res))
    maybe_save(joinpath(prefix, "loss.svg"))

    plot_nbars(params, (r[1][1] for r in res))
    maybe_save(joinpath(prefix, "nbar.svg"))

    plot_hf(params, (r[3][1] for r in res))
    maybe_save(joinpath(prefix, "hf.svg"))
end

const xname = Ref("")
function plot_hook()
    xlabel(xname[])
end

function threadmap(f, arg)
    # A really simple work queue
    n = length(arg)
    counter = Threads.Atomic{Int}(1)
    T = Core.Inference.return_type(f, Tuple{eltype(arg)})
    res = Vector{T}(n)
    nt = Threads.nthreads()
    Threads.@threads for _ in 1:nt
        while true
            i = Threads.atomic_add!(counter, 1)
            if i > n
                break
            end
            res[i] = f(arg[i])
        end
    end
    if !isleaftype(T)
        return [v for v in res]
    end
    return res
end

## End of helper parts

const m_Na = 23f-3 / 6.02f23
const k_Na = Float32(2π) / 589f-9
const k_trap = Float32(2π) / 700f-9
const trap_freq = (67f3, 420f3, 580f3)
η(freq) = Trap.η(m_Na, freq, k_Na)
const η_full = η.(trap_freq)
ηs_Na(a, b, c) = Float32.((a, b, c)) .* η_full
const η_full_trap = Trap.η.(m_Na, trap_freq, k_trap)
ηs_trap(a, b, c) = Float32.((a, b, c)) .* η_full_trap

const η_op = ηs_Na(1, 1, 1)
const η_op_dri = ηs_Na(0, sqrt(0.5), sqrt(0.5))

const δf1 = -25.0e9
const δf2 = -25.0e9 - 1.77e9
const rlof_f1 = (61.542e6 / (δf1 - 1.107266e9))^2
const rlof_f2 = (61.542e6 / (δf2 - 1.107266e9))^2
const rhif_f1 = (61.542e6 / (δf1 + 664.360e6))^2
const rhif_f2 = (61.542e6 / (δf2 + 664.360e6))^2

const rates_f1_coprop = Float32.(all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f1, rlof_f1))
const rates_f1_up = Float32.(all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1))
const rates_f1_down = Float32.(all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1))
const rates_f2_coprop = Float32.(all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f2, rlof_f2))
const rates_f2_counterop = Float32.(all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_f2, rlof_f2))

# Decay rates and Rabi frequencies are measured in MHz (or us⁻¹)
# Times are measured in μs

rates_f1_coprop .*= 4.46e8 / 1e6 # Amp 0.25
rates_f2_coprop .*= 4.1e8 / 1e6 # Amp 0.22
rates_f1_up .*= 1.05e9 / 1e6 # Amp 1.0
rates_f1_down .*= 8.2e8 / 1e6 # Amp 0.22
rates_f2_counterop .*= 3.25e9 / 1e6 # Amp 0.05

function gen_isσ(Δdri)
    res = zeros(Bool, 8, 8)
    idx_to_state = function (idx)
        if idx > 5
            # low F
            return idx - 7
        else
            # high F
            return idx - 3
        end
    end
    for i in 1:8
        mF1 = idx_to_state(i)
        for j in 1:8
            mF2 = idx_to_state(j)
            # From i to j, from mF1 to mF2
            if abs(mF1 + Δdri - mF2) == 1
                res[j, i] = true
            end
        end
    end
    res
end

const isσs_all = [gen_isσ(-1), gen_isσ(0), gen_isσ(1)]

struct BeamSpec{T}
    η::NTuple{3,T}
    rates::Vector{Tuple{Matrix{T},Matrix{Bool}}}
end
function (::Type{BeamSpec{T}})(η, rates::Matrix, pol::NTuple{3,Any}) where T
    rates2 = Tuple{Matrix{T},Matrix{Bool}}[]
    for i in 1:3
        p = pol[i]
        p == 0 && continue
        push!(rates2, (T.(rates .* p), isσs_all[i]))
    end
    return BeamSpec{T}(η, rates2)
end

const bs_f1_coprop = BeamSpec{Float32}(ηs_Na(-sqrt(0.5), 0.5, 0.5), rates_f1_coprop,
                                       (0.5, 0.0, 0.5))
const bs_f2_coprop = BeamSpec{Float32}(ηs_Na(-sqrt(0.5), 0.5, 0.5), rates_f2_coprop,
                                       (0.25, 0.5, 0.25))
const bs_f1_up = BeamSpec{Float32}(ηs_Na(0, sqrt(0.5), -sqrt(0.5)), rates_f1_up,
                                   (0.25, 0.5, 0.25))
const bs_f1_down = BeamSpec{Float32}(ηs_Na(0, -sqrt(0.5), sqrt(0.5)), rates_f1_down,
                                     (0.25, 0.5, 0.25))
const bs_f2_counterop = BeamSpec{Float32}(ηs_Na(0, -sqrt(0.5), -sqrt(0.5)), rates_f2_counterop,
                                          (0.1, 0.0, 0.9))

## Raman powers
# The list of amplitudes we used for Raman transitions
# Co-prop:
#     F1 coprop 0.25 + F2 coprop 0.22
# Axial 1:
#     F1 Up 0.4 + F2 coprop 0.14
#     F1 Up 0.4 + F2 coprop 0.22 (ramp)
#     F1 Up 0.2 + F2 coprop 0.22 (ramp)
#     F1 Up 0.2 + F2 coprop 0.1 (ramp)
# Radial 2:
#     F1 Up 1.0 (ramp) + F2 counterop 0.05
# Radial 3:
#     F1 Down 0.22 (ramp) + F2 counterop 0.05

# Power normalizations:
# We normalize the powers to the ones we calibrated the off-resonance scattering on
# since that is the only single photon measurement we have.
# From there we can calculate the normalized intensities for the conditions listed above.
# Relative powers used for Raman transitions
# Co-prop:
#     F1 coprop 1.000 + F2 coprop 1.000
# Axial 1:
#     F1 Up 0.351 + F2 coprop 0.740
#     F1 Up 0.351 + F2 coprop 1.000 (ramp)
#     F1 Up 0.100 + F2 coprop 1.000 (ramp)
#     F1 Up 0.100 + F2 coprop 0.461 (ramp)
# Radial 2:
#     F1 Up 1.000 (ramp) + F2 counterop 1.000
# Radial 3:
#     F1 Down 1.000 (ramp) + F2 counterop 1.000

# As for the effect of the ramp, since we are doing it only on one beam, the scattering does
# not scale linearly with the instantaneous power.
# Experimentally we see the π time is roughly doubled so we can expect that corresponds to a
# decrease in scattering rate between 2 (square pulse) and 4 (linear ramp).
# Let's just assume the scattering rate is decreased by 3x.

# Now for the calibrated full Rabi frequency with the corresponding normalized powers
# Co-prop:
#     F1 coprop 1.000 + F2 coprop 1.000: Ω = 2π / 26.0us
# Axial 1:
#     F1 Up 0.100 + F2 coprop 0.462: Ω = 2π / 60.2us
# Radial 2:
#     F1 Up 1.000 + F2 counterop 1.000: Ω = 2π / 11.61us
# Radial 3:
#     F1 Down 1.000 + F2 counterop 1.000: Ω = 2π / 11.45us
# The Rabi frequencies (2π times) here are for full matrix element

const sz = 500, 100, 100
const trapscatter = System.Scatter{Float32}(eye(Float32, 8) .* 0.033e-3,
                                            ηs_trap(1, 1, 1), ηs_trap(1, 0, 0),
                                            zeros(Bool, 8, 8),
                                            (0, 1, 1))

function create_raman_raw(t, p1, p2, ramp1, ramp2, bs1::BeamSpec, bs2::BeamSpec, Ω0, Δn)
    Ω = sqrt(p1 * p2) * Ω0
    if ramp1
        p1 /= 3
        Ω /= 2
        @assert !ramp2
    end
    if ramp2
        p2 /= 3
        Ω /= 2
    end
    s1s = [System.Scatter{Float32}(p1 * r[1], η_op, abs.(bs1.η), r[2],
                                   (0, 1, 1)) for r in bs1.rates]
    s2s = [System.Scatter{Float32}(p2 * r[1], η_op, abs.(bs2.η), r[2],
                                   (0, 1, 1)) for r in bs2.rates]
    return System.RealRaman{Float32,1,6}(t, Ω, abs.(bs1.η .- bs2.η), Δn, sz,
                                         [s1s; s2s; [trapscatter]])
end

struct RamanSpec{T}
    Ω0::T
    bs1::BeamSpec{T}
    bs2::BeamSpec{T}
end

const raman_specs = [RamanSpec{Float32}(0.2417, bs_f1_coprop, bs_f2_coprop), # co-prop
                     RamanSpec{Float32}(0.4856, bs_f1_up, bs_f2_coprop), # axial 1
                     RamanSpec{Float32}(0.5412, bs_f1_up, bs_f2_counterop), # radial 2
                     RamanSpec{Float32}(0.5487, bs_f1_down, bs_f2_counterop) # radial 3
                     ]

get_Δns(ax, Δn::NTuple{3,Any}) = Δn
function get_Δns(ax, Δn)
    if ax == 0
        @assert Δn == 0
        return (0, 0, 0)
    elseif ax == 1
        return (Δn, 0, 0)
    elseif ax == 2
        return (0, Δn, 0)
    elseif ax == 3
        return (0, 0, Δn)
    else
        error("Invalid axis number $ax")
    end
end

function create_raman(t, p1, p2, ramp1, ramp2, ax, Δn=0)
    rs = raman_specs[ax + 1]
    return create_raman_raw(t, p1, p2, ramp1, ramp2, rs.bs1, rs.bs2, rs.Ω0, get_Δns(ax, Δn))
end
create_wait(t) = System.MultiOP{Float32}(t, [trapscatter])

const BuilderT = Setup.SeqBuilder{System.StateC,Void}
statec() = System.StateC(sz...)

function create_sequence(t)
    builder = BuilderT(System.ThermalInit{1,Float32}(0, 0, 0), Setup.Dummy(),
                       Setup.CombinedMeasure(System.NBarMeasure(),
                                             System.GroundStateMeasure(),
                                             System.HyperFineMeasure{8}()))
    Setup.add_pulse(builder, create_wait(t * 1e3))
    Setup.add_pulse(builder, create_raman(81, 0.100, 0.461, false, true, 1, -1))
    return builder.seq
end

const params = linspace(0, 300, 100)
res = @time threadmap(p->Setup.run(create_sequence(p), statec(), nothing, 100000), params)

if interactive()
    plot_result(params, res, "")
else
    mkpath(ARGS[1], 0o750)
    plot_result(params, res, ARGS[1])
end
maybe_show()
