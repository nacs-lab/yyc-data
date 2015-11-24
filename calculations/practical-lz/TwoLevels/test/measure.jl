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
    y = Complex64[1, 2]
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
    measures = Any[(1, 10_000)=>DummyMeasure,
                   (20_000, 30_000)=>DummyMeasure,
                   (50_000, 3_000_000)=>DummyMeasure,
                   (5_000_000, 90_000_000)=>DummyMeasure]
    dt = 1f-3
    n = 100_000_000
    measure_list = MeasureList{Float32}(1:n, dt, measures)
    y = Complex64[1, 2]

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
    y = Complex64[1, 2]
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
    imin::Int
    imax::Int
    idx::Int
    CountMeasure(idxs, dt, idx=1) = new(0, first(idxs), last(idxs), idx)
end
function Measures.snapshot(c::CountMeasure, y, idx, t)
    c.n += 1
    @assert idx == c.n
end
function Base.reset(c::CountMeasure)
    c.n = 0
end
function test_list()
    measures = Any[(1, 10_000)=>CountMeasure,
                   (20_000, 30_000)=>(CountMeasure, 2),
                   (50_000, 3_000_000)=>(CountMeasure, 3),
                   (5_000_000, 90_000_000)=>(CountMeasure, 4)]
    dt = 1f-3
    n = 100_000_000
    measure_list = MeasureList{Float32}(1:n, dt, measures)
    y = Complex64[1, 2]

    # @code_llvm Measures.snapshot(measure_list, y, 0, 0 * dt)
    Measures.snapshot(measure_list, y, 1, 1 * dt)
    t = @elapsed for i in 2:(n + 1)
        Measures.snapshot(measure_list, y, i, i * dt)
    end
    info(@sprintf("    Time per measure: %.2fns", t / n * 1e9))
    for (i, ((imin, imax), m)) in enumerate(measure_list)
        measure = m.measure::CountMeasure{Float32}
        @test i == measure.idx
        @test imin == measure.imin
        @test imax == measure.imax
        @test imax - imin + 1 == measure.n
    end
    reset(measure_list)
    for ((imin, imax), m) in measure_list
        measure = m.measure::CountMeasure{Float32}
        @test measure.n == 0
    end
end
test_list()

end
