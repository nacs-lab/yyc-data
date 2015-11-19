#!/usr/bin/julia -f

module TestMeasure

using Base.Test
import TwoLevels: AbstractMeasure, Measures
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
    #     Measures.snapshot(dummy_measure, y, i, i * dt)
    # end
    # @code_llvm Measures.snapshot(wrapped_measure, y, 0, 0 * dt)
    # @code_native Measures.snapshot(wrapped_measure, y, 0, 0 * dt)
    Measures.snapshot(wrapped_measure, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        Measures.snapshot(wrapped_measure, y, i, i * dt)
    end
    info(@sprintf("    Time per measure: %.2fns", t / n * 1e9))
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

    # @code_llvm Measures.snapshot(measure_list, y, 0, 0 * dt)
    Measures.snapshot(measure_list, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        Measures.snapshot(measure_list, y, i, i * dt)
    end
    info(@sprintf("    Time per measure: %.2fns", t / n * 1e9))
end
perf_list()

info("  Test FullMeasure performance")
function perf_full()
    n = 10_000_000

    full_measure = FullMeasure{Float32}(n)
    wrapped_measure = MeasureWrapper(full_measure)
    y = Float32[1, 2]
    dt = 1f-3

    Measures.snapshot(full_measure, y, 1, 1 * dt)
    # @code_llvm Measures.snapshot(full_measure, y, 0, 0 * dt)
    t = @elapsed for i in 2:(n + 1)
        Measures.snapshot(full_measure, y, i, i * dt)
    end
    info(@sprintf("    Time per full measure: %.2fns", t / n * 1e9))
    Measures.snapshot(wrapped_measure, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        Measures.snapshot(wrapped_measure, y, i, i * dt)
    end
    info(@sprintf("    Time per wrapped full measure: %.2fns", t / n * 1e9))
end
perf_full()

info("  Test MeasureList")
type CountMeasure{T} <: AbstractMeasure{T}
    n::Int
    CountMeasure() = new(0)
end
function Measures.snapshot(c::CountMeasure, y, idx, t)
    c.n += 1
    @assert idx == c.n
end
function Measures.reset(c::CountMeasure)
    c.n = 0
end
function test_list()
    TM = Pair{Tuple{Int,Int},MeasureWrapper{Float32}}
    measures = TM[(1, 10_000)=>MeasureWrapper(CountMeasure{Float32}()),
                  (20_000, 30_000)=>MeasureWrapper(CountMeasure{Float32}()),
                  (50_000, 3_000_000)=>MeasureWrapper(CountMeasure{Float32}()),
                  (5_000_000, 90_000_000)=>MeasureWrapper(CountMeasure{Float32}())]
    dt = 1f-3
    measure_list = MeasureList{Float32}(measures, dt)
    n = 100_000_000
    y = Float32[1, 2]

    # @code_llvm Measures.snapshot(measure_list, y, 0, 0 * dt)
    Measures.snapshot(measure_list, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        Measures.snapshot(measure_list, y, i, i * dt)
    end
    info(@sprintf("    Time per measure: %.2fns", t / n * 1e9))
    for ((imin, imax), m) in measures
        @test imax - imin + 1 == m.measure.n
    end
    Measures.reset(measure_list)
    for ((imin, imax), m) in measures
        @test m.measure.n == 0
    end
end
test_list()

end
