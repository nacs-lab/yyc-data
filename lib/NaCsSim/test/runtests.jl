#!/usr/bin/julia

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
state = System.StateC(sz...)
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

function add_group1(builder)
    add_pulse(builder,
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 6, 8),
                          RamanParams(1, 5, 10),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          12))
end
function add_group2(builder)
    add_pulse(builder,
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 5, 10),
                          RamanParams(1, 4, 10),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          12))
end
function add_group3(builder)
    add_pulse(builder,
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 4, 10),
                          RamanParams(1, 3, 12),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          12))
end
function add_group4(builder)
    add_pulse(builder,
              Grp2AParams(OPParams(15, 0.3, 0.01),
                          RamanParams(1, 3, 12),
                          RamanParams(1, 2, 4),
                          RamanParams(2, 1, 5),
                          RamanParams(3, 1, 5),
                          12))
end
function add_group5(builder)
    add_pulse(builder,
              Grp2AParams(OPParams(15, 0.06, 0.01),
                          RamanParams(1, 2, 4),
                          RamanParams(1, 1, 4),
                          RamanParams(2, 1, 3),
                          RamanParams(3, 1, 3),
                          50))
end

# builder = BuilderT(init, Setup.Dummy(), System.HyperFineMeasure{3}())
# builder = BuilderT(init, Setup.Dummy(), System.NBarMeasure())
# builder = BuilderT(init, Setup.Dummy(), System.GroundStateMeasure())
# builder = BuilderT(init, Setup.Dummy(), Setup.Dummy())
# add_group1(builder)
# add_group2(builder)
# add_group3(builder)
# add_group4(builder)
# add_group5(builder)

function create_sequence(ngroup, nbar)
    if nbar
        builder = BuilderT(init, Setup.Dummy(), System.NBarMeasure())
    else
        builder = BuilderT(init, Setup.Dummy(), System.GroundStateMeasure())
    end
    # builder = BuilderT(init, Setup.Dummy(), Setup.Dummy())
    # builder = BuilderT(init, Setup.Dummy(), System.HyperFineMeasure{3}())
    adders = [add_group1, add_group2, add_group3, add_group4, add_group5]
    for i in 1:ngroup
        adders[i](builder)
    end
    return builder.seq
end

function run_sequences()
    for i in 0:5
        println("Pulse groups: $i")
        nbars, = Setup.run(create_sequence(i, true), state, nothing, 100000)
        ground = Setup.run(create_sequence(i, false), state, nothing, 100000)
        println("    nbar: $nbars; pgrd: $ground")
    end
end

# run_sequences()

seq1 = create_sequence(5, true)
seq2 = create_sequence(5, false)

# @time @show Setup.run(seq1, state, nothing, 10)
# @time @show Setup.run(seq2, state, nothing, 10)

@time @show Setup.run(seq1, state, nothing, 100000)
@time @show Setup.run(seq2, state, nothing, 100000)
