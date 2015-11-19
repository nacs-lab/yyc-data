#!/usr/bin/julia -f

module Measures

abstract AbstractMeasure{T}

function measure_snapshot end

immutable DummyMeasure{T} <: AbstractMeasure{T}
end

@inline measure_snapshot(::DummyMeasure, y, idx, t) = nothing

function _measure_wrapper{T}(measure::AbstractMeasure{T}, y::Vector{T},
                             idx::Int, t::T)
    measure_snapshot(measure, y, idx, t)
    nothing
end

# Efficient type stable wrapper of arbitrary AbstractMeasure
# The overhead of calling a wrapped AbstractMeasure is roughly two function
# calls and a GC frame push and pop (plus no inlining). This is measured to be
# ~4-5ns on my laptop.
immutable MeasureWrapper{T} <: AbstractMeasure{T}
    fptr::Ptr{Void}
    measure::AbstractMeasure{T} # For GC root
    function MeasureWrapper(measure::AbstractMeasure{T})
        M = typeof(measure)
        fptr = cfunction(_measure_wrapper, Void,
                         Tuple{Ref{M},Ref{Vector{T}},Int,T})
        new(fptr, measure)
    end
end
MeasureWrapper{T}(measure::AbstractMeasure{T}) = MeasureWrapper{T}(measure)

@inline function measure_snapshot{T}(wrapper::MeasureWrapper{T}, y::Vector{T},
                                     idx::Int, t::T)
    # This let LLVM optimize out the jl_throw when fptr is NULL
    # Should be replaced with `llvm.assume` once I upgrade my julia version
    wrapper.fptr == C_NULL && return
    ccall(wrapper.fptr, Void, (Any, Any, Int, T), wrapper.measure, y, idx, t)
end

type MeasureList{T} <: AbstractMeasure{T}
    cur_measure::MeasureWrapper{T}
    tidx_max::Int
    tidx_offset::Int
    t_offset::T
    next_measure::Int
    # Following fields should not be mutated
    measures::Vector{Pair{Tuple{Int,Int},MeasureWrapper{T}}}
    dummy_measure::MeasureWrapper{T}
    dt::T
    function MeasureList(measures, dt)
        # TODO verify if the measures list is valid
        dummy_measure = MeasureWrapper{T}(DummyMeasure{T}())
        new(dummy_measure, 0, 0, zero(T), isempty(measures) ? 0 : 1,
            measures, dummy_measure, dt)
    end
end

@noinline function measure_list_update(list, idx)
    list.tidx_offset = idx - 1
    list.t_offset = (idx - 1) * list.dt
    next_midx = list.next_measure
    if next_midx == 0
        list.cur_measure = list.dummy_measure
        list.tidx_max = typemax(Int)
        return
    end
    @inbounds _next = list.measures[next_midx]
    next_tidx_min, next_tidx_max = _next.first
    if idx < next_tidx_min
        list.cur_measure = list.dummy_measure
        list.tidx_max = next_tidx_min - 1
    else
        list.cur_measure = _next.second
        list.tidx_max = next_tidx_max
        list.next_measure = ifelse(next_midx >= length(list.measures),
                                   0, next_midx + 1)
    end
    nothing
end

@inline function measure_snapshot{T}(list::MeasureList{T}, y::Vector{T},
                                     idx::Int, t::T)
    list.tidx_max < idx && measure_list_update(list, idx)
    measure_snapshot(list.cur_measure, y, idx - list.tidx_offset,
                     t - list.t_offset)
end

immutable FullMeasure{T} <: AbstractMeasure{T}
    ys::Matrix{T}
    FullMeasure(nsteps) = new(Matrix{T}(2, nsteps + 1))
end

@inline function measure_snapshot(measure::FullMeasure, y, idx, t)
    @inbounds measure.ys[1, idx] = abs2(y[1])
    @inbounds measure.ys[2, idx] = abs2(y[2])
    nothing
end

end
