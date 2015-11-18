#!/usr/bin/julia -f

module TestMeasure

import TwoLevels: AbstractMeasure, measure_snapshot
import TwoLevels: DummyMeasure, MeasureWrapper

info("Test Measure")

info("  Test MeasureWrapper performance")
function perf()
    dummy_measure = DummyMeasure()
    wrapped_measure = MeasureWrapper{Float32}(dummy_measure)
    n = 100_000_000
    y = Float32[1, 2]
    dt = 1f-3

    # @time for i in 1:n
    #     measure_snapshot(dummy_measure, y, i, i * dt)
    # end
    # @code_llvm measure_snapshot(wrapped_measure, y, 0, 0 * dt)
    measure_snapshot(wrapped_measure, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        measure_snapshot(wrapped_measure, y, i, i * dt)
    end
    info(@sprintf("    Time per wrapper run: %.1fns", t / n * 1e9))
end

perf()

end
