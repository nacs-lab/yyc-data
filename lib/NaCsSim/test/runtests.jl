#!/usr/bin/julia

import NaCsSim: Setup, System
import NaCsCalc: Trap

const BuilderT = Setup.SeqBuilder{System.State{Float32,3},Void}
const sz = 500, 100, 100
const init = System.ThermalInit{1,Float32}(20, 4, 4)

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

# 1: (2, -2); 2: (2, -1); 3: (1, -1)
# builder = BuilderT(init, Setup.Dummy(), System.HyperFineMeasure())
builder = BuilderT(init, Setup.Dummy(), System.NBarMeasure())
# builder = BuilderT(init, Setup.Dummy(), Setup.Dummy())
state = System.State{Float32,3}(sz...)
function f1op(t, γ, ηs=η_op, ηdri=η_op_dri)
    System.OP{Float32}(t, branching_11 .* γ, ηs, ηdri, isσ)
end
function f2op(t, γ1, γ2, γ2′=0f0 * γ2, ηs=η_op, ηdri=η_op_dri)
    branching = (branching_11 .* γ1 .+ branching_21 .* γ2 .+
                 branching_22 .* γ2′)
    System.OP{Float32}(t, branching, ηs, ηdri, isσ)
end
raman_pulse1_6 = System.Raman{Float32,1,3}(8, 1, η_raman1, (-6, 0, 0), sz)
raman_pulse2_2 = System.Raman{Float32,1,3}(5, 1, η_raman2, (0, -2, 0), sz)
raman_pulse3_2 = System.Raman{Float32,1,3}(5, 1, η_raman3, (0, 0, -2), sz)
f2op_pulse = f2op(40, 5, 0.5)

for i in 1:100
    Setup.add_pulse(builder, raman_pulse1_6)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse2_2)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse1_6)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse3_2)
    Setup.add_pulse(builder, f2op_pulse)
end
@time Setup.run(builder.seq, state, nothing, 1)
@show @time Setup.run(builder.seq, state, nothing, 10000)
