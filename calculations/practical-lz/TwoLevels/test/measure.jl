#!/usr/bin/julia -f

module TestMeasure

import TwoLevels: AbstractMeasure, measure_snapshot
import TwoLevels: MeasureWrapper, MeasureList
import TwoLevels: DummyMeasure, FullMeasure

info("Test Measure")

info("  Test MeasureWrapper performance")
function perf_wrapper()
    dummy_measure = DummyMeasure{Float32}()
    wrapped_measure = MeasureWrapper(dummy_measure)
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
    info(@sprintf("    Time per measure: %.1fns", t / n * 1e9))
end
perf_wrapper()

info("  Test MeasureList performance")
function perf_list()
    dummy_measure = DummyMeasure{Float32}()
    TM = Pair{Tuple{Int,Int},MeasureWrapper{Float32}}
    measures = TM[(1, 10_000)=>MeasureWrapper(dummy_measure),
                  (20_000, 30_000)=>MeasureWrapper(dummy_measure),
                  (50_000, 3_000_000)=>MeasureWrapper(dummy_measure),
                  (5_000_000, 90_000_000)=>MeasureWrapper(dummy_measure)]
    dt = 1f-3
    measure_list = MeasureList{Float32}(measures, dt)
    n = 100_000_000
    y = Float32[1, 2]

    # @code_llvm measure_snapshot(measure_list, y, 0, 0 * dt)
    measure_snapshot(measure_list, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        measure_snapshot(measure_list, y, i, i * dt)
    end
    info(@sprintf("    Time per measure: %.1fns", t / n * 1e9))
end
perf_list()

info("  Test FullMeasure performance")
function perf_full()
    n = 100_000_000

    full_measure = FullMeasure{Float32}(100_000_000)
    wrapped_measure = MeasureWrapper(full_measure)
    y = Float32[1, 2]
    dt = 1f-3

    measure_snapshot(full_measure, y, 1, 1 * dt)
    # @code_llvm measure_snapshot(full_measure, y, 0, 0 * dt)
    t = @elapsed for i in 2:(n + 1)
        measure_snapshot(full_measure, y, i, i * dt)
    end
    info(@sprintf("    Time per full measure: %.1fns", t / n * 1e9))
    measure_snapshot(wrapped_measure, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        measure_snapshot(wrapped_measure, y, i, i * dt)
    end
    info(@sprintf("    Time per wrapped full measure: %.1fns", t / n * 1e9))
end
perf_full()

end
