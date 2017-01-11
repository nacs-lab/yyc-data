#!/usr/bin/julia

import NaCsSim: Setup, System

const BuilderT = Setup.SeqBuilder{System.State{Float32,3},Void}

const sz = 500, 100, 100

const init = System.ThermalInit{1,Float32}(20, 4, 4)

# 1: (2, -2); 2: (2, -1); 3: (1, -1)
# builder = BuilderT(init, Setup.Dummy(), System.HyperFineMeasure())
builder = BuilderT(init, Setup.Dummy(), System.NBarMeasure())
# builder = BuilderT(init, Setup.Dummy(), Setup.Dummy())
state = System.State{Float32,3}(sz...)
function f1op(γ, t, ηs=(0.65f0, 0.24f0, 0.24f0), ηdri=(0f0, 0.17f0, 0.17f0))
    ratio = Float32[0 0 1 / 6
                    0 0 1 / 12
                    0 0 1 / 4]
    isσ = [false false false
            false false true
            false false true]
    System.OP{Float32}(t, ratio * γ, ηs, ηdri, isσ)
end
function f2op(γ1, γ2, t, ηs=(0.0065f0, 0.0024f0, 0.0024f0), ηdri=(0f0, 0.0017f0, 0.0017f0))
    ratio1 = Float32[0 0 1 / 6
                     0 0 1 / 12
                     0 0 1 / 4]
    ratio2 = Float32[0 1 / 3 0
                     0 1 / 6 0
                     0 1 / 2 0]
    ratio2′ = Float32[1 / 3 0 0
                       1 / 6 0 0
                       1 / 2 0 0] * 0.01f0
    isσ = [false false false
            true true true
            true true true]
    System.OP{Float32}(t, ratio1 * γ1 + (ratio2 + ratio2′) * γ2, ηs, ηdri, isσ)
end
raman_pulse = System.Raman{Float32,1,3}(5, 1, (0.0f0, 0.35f0, 0f0),
                                        (0, -2, 0), sz)
raman_pulse2 = System.Raman{Float32,1,3}(5, 1, (0.0f0, 0f0, 0.35f0),
                                         (0, 0, -2), sz)
raman_pulse3 = System.Raman{Float32,1,3}(10, 1, (0.45f0, 0f0, 0f0),
                                         (-6, 0, 0), sz)
f2op_pulse = f2op(5, 0.5, 40)

for i in 1:100
    Setup.add_pulse(builder, raman_pulse3)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse3)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse2)
    Setup.add_pulse(builder, f2op_pulse)
    Setup.add_pulse(builder, raman_pulse3)
    Setup.add_pulse(builder, f2op_pulse)
end
@time Setup.run(builder.seq, state, nothing, 1)
@show @time Setup.run(builder.seq, state, nothing, 10000)
