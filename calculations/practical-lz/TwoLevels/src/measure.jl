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

end
