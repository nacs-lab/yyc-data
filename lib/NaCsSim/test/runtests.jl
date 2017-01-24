#!/usr/bin/julia

@everywhere import NaCsSim: Setup, System

@everywhere module TestSequence
import NaCsSim: Setup, System
import NaCsCalc: Trap

const BuilderT = Setup.SeqBuilder{System.StateC,Void}
const sz = 500, 100, 100
const init = System.ThermalInit{1,Float32}(15, 4, 4)

const m_Na = 23f-3 / 6.02f23
const k_Na = Float32(2π) / 589f-9
η(freq) = Trap.η(m_Na, freq, k_Na)

const η_op = (η(60f3), η(420f3), η(420f3))
const η_op_dri = (0f0, η(420f3) / sqrt(2f0), η(420f3) / sqrt(2f0))
const isσ = [false false false
              true true true
              true true true]
const branching_11 = Float32[0 0 1 / 6
                             0 0 1 / 12
                             0 0 1 / 4]
const branching_21 = Float32[0 1 / 3 0
                             0 1 / 6 0
                             0 1 / 2 0]
const branching_22 = Float32[1 / 3 0 0
                             1 / 6 0 0
                             1 / 2 0 0]
const η_raman1 = (η(60f3) / sqrt(2f0), 0f0, 0f0)
const η_raman2 = (0f0, η(420f3) * sqrt(2f0), 0f0)
const η_raman3 = (0f0, 0f0, η(420f3) * sqrt(2f0))
const η_ramans = (η_raman1, η_raman2, η_raman3)

# 1: (2, -2); 2: (2, -1); 3: (1, -1)
export statec
const statec = System.StateC(sz...)
function op_pulse(t, γ1, γ2, γ2′=0.01f0 * γ2, ηs=η_op, ηdri=η_op_dri)
    branching = (branching_11 .* γ1 .+ branching_21 .* γ2 .+
                 branching_22 .* γ2′)
    System.OP{Float32}(t, branching, ηs, ηdri, isσ)
end
function raman_pulse(ax, order, t, Ω=1)
    ns = ntuple(i->i == ax ? -order : 0, Val{3})
    System.Raman{Float32,1,3}(t, Ω, η_ramans[ax], ns, sz)
end

immutable OPParams
    γ1::Float32
    γ2::Float32
    darkness::Float32
end

pulse(params::OPParams) =
    op_pulse(1, params.γ1, params.γ2, params.γ2 * params.darkness)

immutable RamanParams
    ax::Int
    order::Int
    t::Float32
end

pulse(params::RamanParams) = raman_pulse(params.ax, params.order, params.t)

add_pulse(builder, params) = Setup.add_pulse(builder, pulse(params))

# Group with two axial cooling per loop
immutable Grp2AParams
    op::OPParams
    raman11::RamanParams
    raman12::RamanParams
    raman2::RamanParams
    raman3::RamanParams
    rep::Int
end

function add_pulse(builder, params::Grp2AParams)
    op = pulse(params.op)
    for i in 1:params.rep
        Setup.add_pulse(builder, pulse(params.raman11))
        Setup.add_pulse(builder, op)
        Setup.add_pulse(builder, pulse(params.raman12))
        Setup.add_pulse(builder, op)
        Setup.add_pulse(builder, pulse(params.raman2))
        Setup.add_pulse(builder, op)
        Setup.add_pulse(builder, pulse(params.raman11))
        Setup.add_pulse(builder, op)
        Setup.add_pulse(builder, pulse(params.raman12))
        Setup.add_pulse(builder, op)
        Setup.add_pulse(builder, pulse(params.raman3))
        Setup.add_pulse(builder, op)
    end
end

export create_sequence
function create_sequence(ngroup, ncycles)
    pulses_left = Ref(ncycles)

    take_pulses = function (n)
        nleft = pulses_left[]
        if nleft > n
            nleft -= n
            pulses_left[] = nleft
            return n
        else
            pulses_left[] = 0
            return nleft
        end
    end

    # builder = BuilderT(init, Setup.Dummy(), Setup.Dummy())
    # builder = BuilderT(init, Setup.Dummy(), System.HyperFineMeasure{3}())
    builder = BuilderT(init, Setup.Dummy(),
                       Setup.CombinedMeasure(System.NBarMeasure(),
                                             System.GroundStateMeasure()))

    pulses = [Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 6, 8),
                          RamanParams(1, 5, 10),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 5, 10),
                          RamanParams(1, 4, 10),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 4, 10),
                          RamanParams(1, 3, 12),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 3, 12),
                          RamanParams(1, 2, 4),
                          RamanParams(2, 1, 5),
                          RamanParams(3, 1, 5),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.06, 0.01),
                          RamanParams(1, 2, 4),
                          RamanParams(1, 1, 4),
                          RamanParams(2, 1, 3),
                          RamanParams(3, 1, 3),
                          take_pulses(50))]
    for i in 1:ngroup
        add_pulse(builder, pulses[i])
    end
    return builder.seq
end

end

@everywhere using TestSequence

const params = 0:98

res = pmap(ncycles->Setup.run(create_sequence(5, ncycles), statec,
                              nothing, 100000), params)

using PyPlot

function plot_result(params, res)
    figure()
    gp = [r[2].a for r in res]
    gp_unc = [r[2].s for r in res]
    errorbar(params, gp, gp_unc, label="Ground state")

    nbar_res = (r[1] for r in res)
    total_res = (nr[2] for nr in nbar_res)
    total = [t.a for t in total_res]
    total_unc = [t.s for t in total_res]
    errorbar(params, total, total_unc, label="Total")
    legend()
    grid()

    figure()
    nbarx_res = (nr[1][1] for nr in nbar_res)
    nbary_res = (nr[1][2] for nr in nbar_res)
    nbarz_res = (nr[1][3] for nr in nbar_res)

    nbarx = [n.a for n in nbarx_res]
    nbarx_unc = [n.s for n in nbarx_res]
    nbary = [n.a for n in nbary_res]
    nbary_unc = [n.s for n in nbary_res]
    nbarz = [n.a for n in nbarz_res]
    nbarz_unc = [n.s for n in nbarz_res]
    gca()[:set_yscale]("log")
    errorbar(params, nbarx, nbarx_unc, label="X")
    errorbar(params, nbary, nbary_unc, label="Y")
    errorbar(params, nbarz, nbarz_unc, label="Z")
    legend()
    grid()
end

plot_result(params, res)
show()

# function run_sequences()
#     for i in 0:5
#         println("Pulse groups: $i")
#         (nbars,), ground =
#             Setup.run(create_sequence(i), statec, nothing, 100000)
#         println("    nbar: $nbars; pgrd: $ground")
#     end
# end
