#!/usr/bin/julia -f

module TestMeasure

import TwoLevels: AbstractMeasure, measure_snapshot
import TwoLevels: DummyMeasure, MeasureWrapper

function perf()
    dummy_measure = DummyMeasure()
    wrapped_measure = MeasureWrapper{Float32}(dummy_measure)
    n = 1000_000
    y = Float32[1, 2]
    dt = 1f-3

    @time for i in 1:n
        measure_snapshot(dummy_measure, y, i, i * dt)
    end
    # @code_llvm measure_snapshot(wrapped_measure, y, 0, 0 * dt)
    @time for i in 1:n
        measure_snapshot(wrapped_measure, y, i, i * dt)
    end
end

perf()

end
