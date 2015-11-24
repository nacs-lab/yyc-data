#!/usr/bin/julia -f

module TestSequence

import TwoLevels: Sequence, propagate
import TwoLevels: DummyMeasure, MeasureList
import TwoLevels: ConstDrive

info("  Test propagate(::Sequence) performance")
function perf_seq()
    dummy_measure = DummyMeasure{Float32}()
    dt = 1f-3
    n = 100_000_000
    y = Float32[1, 0]

    seq = Sequence{Float32}(dt, n, (n=>ConstDrive(1f0, 2f0),), dummy_measure)
    t = @elapsed propagate(seq, y)
    info(@sprintf("    Time per iteration: %.2fns", t / n * 1e9))
end
perf_seq()

end
