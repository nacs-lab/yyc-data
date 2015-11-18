#!/usr/bin/julia -f

module Measure

abstract AbstractMeasure

function measure_snapshot end

immutable DummyMeasure <: AbstractMeasure
end

@inline function measure_snapshot(::DummyMeasure, y, idx, t)
    nothing
end

function _measure_wrapper{T}(measure, y::Vector{T}, idx::Int, t::T)
    measure_snapshot(measure, y, idx, t)
    nothing
end

# Efficient type stable wrapper of arbitrary AbstractMeasure
# The overhead of calling a wrapped AbstractMeasure is roughly two function
# calls and a GC frame push and pop (plus no inlining). This is measured to be
# ~4-5ns on my laptop.
immutable MeasureWrapper{T} <: AbstractMeasure
    fptr::Ptr{Void}
    pmeasure::Ptr{Void} # To avoid the UndefRef check
    measure # For GC root
    function MeasureWrapper{M}(measure::M)
        fptr = cfunction(_measure_wrapper, Void,
                         Tuple{Ref{M},Ref{Vector{T}},Int,T})
        # This is not compatible with moving GC
        new(fptr, pointer_from_objref(measure), measure)
    end
end

@inline function measure_snapshot{T}(wrapper::MeasureWrapper{T}, y::Vector{T},
                                     idx::Int, t::T)
    # This let LLVM optimize out the jl_throw when fptr is NULL
    # Should be replaced with `llvm.assume` once I upgrade my julia version
    wrapper.fptr == C_NULL && return
    ccall(wrapper.fptr, Void, (Ptr{Void}, Any, Int, T),
          wrapper.pmeasure, y, idx, t)
end

type MeasureList{T} <: AbstractMeasure
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
        dummy_measure = MeasureWrapper{T}(DummyMeasure())
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

end
