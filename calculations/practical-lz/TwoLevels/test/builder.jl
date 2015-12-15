#!/usr/bin/julia -f

module TestBuilder

import TwoLevels: SeqBuilder, Builders
import TwoLevels: propagate
import TwoLevels: DummyMeasure, MeasureList, FullMeasure
import TwoLevels: ConstDrive

info("  Test seq builder")
function perf_builder()
    dt = 1f-3
    n = 100_000_000
    builder = SeqBuilder(dt)

    # Builders.start_measure(builder, FullMeasure)
    Builders.start_measure(builder, DummyMeasure)
    Builders.add_drive(builder, n * dt, ConstDrive(1f0, 2f0))
    Builders.finish_measure(builder)

    y = Float32[1, 0]
    seq = Builders.build(builder)
    t = @elapsed propagate(seq, y)
    info(@sprintf("    Time per iteration: %.2fns", t / n * 1e9))
end
perf_builder()

end
