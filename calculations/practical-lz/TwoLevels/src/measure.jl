#!/usr/bin/julia -f

module Measures

if VERSION >= v"0.5.0-dev+1297"
    @inline assume(x::Bool) =
        Base.llvmcall(("declare void @llvm.assume(i1)",
                       """call void @llvm.assume(i1 %0)
                       ret void"""), Void, Tuple{Bool}, x)
else
    @inline assume(x::Bool) = nothing
end

abstract AbstractMeasure{T}

function snapshot end
Base.reset(::AbstractMeasure) = nothing

immutable DummyMeasure{T} <: AbstractMeasure{T}
    DummyMeasure() = new()
    # Standard builder interface
    DummyMeasure(idxs, dt) = DummyMeasure{T}()
end

@inline snapshot(::DummyMeasure, y, idx, t) = nothing

function _measure_wrapper{T}(measure::AbstractMeasure{T},
                             y::Vector{Complex{T}}, idx::Int, t::T)
    snapshot(measure, y, idx, t)
    nothing
end

# Efficient type stable wrapper of arbitrary AbstractMeasure
# The overhead of calling a wrapped AbstractMeasure is roughly two function
# calls and a GC frame push and pop (plus no inlining). This is measured to be
# ~4-5ns on my laptop.
immutable MeasureWrapper{T} <: AbstractMeasure{T}
    fptr::Ptr{Void}
    measure
    function MeasureWrapper(measure::AbstractMeasure{T})
        M = typeof(measure)
        fptr = cfunction(_measure_wrapper, Void,
                         Tuple{Ref{M},Ref{Vector{Complex{T}}},Int,T})
        new(fptr, measure)
    end
    # Standard builder interface
    MeasureWrapper(idxs, dt, measure::AbstractMeasure{T}) =
        MeasureWrapper{T}(measure)
    MeasureWrapper{M<:AbstractMeasure}(idxs, dt, ::Type{M}, args...) =
        MeasureWrapper{T}(M{T}(idxs, dt, args...))
end
MeasureWrapper{T}(measure::AbstractMeasure{T}) = MeasureWrapper{T}(measure)
Base.reset(wrapper::MeasureWrapper) = reset(wrapper.measure)

@inline function snapshot{T}(wrapper::MeasureWrapper{T}, y::Vector{Complex{T}},
                             idx::Int, t::T)
    fptr = wrapper.fptr
    assume(fptr != C_NULL)
    ccall(fptr, Void, (Any, Any, Int, T), wrapper.measure, y, idx, t)
end

typealias MeasuresVector{T} Vector{Pair{Tuple{Int,Int},MeasureWrapper{T}}}

function verify_measure_list{T}(measures::MeasuresVector{T})
    prev_imax = 0
    @inbounds for ((imin, imax), m) in measures
        imin <= prev_imax && ArgumentError("Overlapping measures")
        imax < imin && ArgumentError("Measure of negative length")
    end
end

function to_measure_list_item{T}(::Type{T}, _idxs, imin, imax, dt, m::ANY)
    idxs = _idxs[imin]:_idxs[imax]
    if isa(m, MeasureWrapper{T})
        return m::MeasureWrapper{T}
    elseif isa(m, AbstractMeasure{T})
        return MeasureWrapper{T}(m)::MeasureWrapper{T}
    elseif isa(m, DataType) && (m::DataType) <: AbstractMeasure
        return MeasureWrapper{T}(idxs, dt, m)::MeasureWrapper{T}
    end
    return MeasureWrapper{T}(idxs, dt, m...)::MeasureWrapper{T}
end

type MeasureList{T} <: AbstractMeasure{T}
    cur_measure::MeasureWrapper{T}
    tidx_max::Int
    tidx_offset::Int
    t_offset::T
    next_measure::Int
    # Following fields should not be mutated
    measures::MeasuresVector{T}
    dummy_measure::MeasureWrapper{T}
    dt::T
    function MeasureList(measures::MeasuresVector{T}, dt)
        verify_measure_list(measures)
        dummy_measure = MeasureWrapper{T}(DummyMeasure{T}())
        new(dummy_measure, 0, 0, zero(T), isempty(measures) ? 0 : 1,
            measures, dummy_measure, dt)
    end
    MeasureList(idxs, dt, measures) =
        MeasureList{T}([Pair((imin::Int, imax::Int),
                             to_measure_list_item(T, idxs, imin, imax,
                                                  dt, m)::MeasureWrapper{T})
                        for ((imin, imax), m) in measures], dt)
end
function Base.reset{T}(list::MeasureList{T})
    list.cur_measure = list.dummy_measure
    list.tidx_max = 0
    list.tidx_offset = 0
    list.t_offset = zero(T)
    list.next_measure = isempty(list.measures) ? 0 : 1
    @inbounds for ((imin, imax), m) in list.measures
        reset(m)
    end
end
function Base.getindex(list::MeasureList, idx)
    idx, m = list.measures[idx]
    idx=>(m.measure::AbstractMeasure)
end
Base.length(list::MeasureList) = length(list.measures)
Base.start(list::MeasureList) = 1
Base.next(list::MeasureList, idx) = (list[idx], idx + 1)
Base.done(list::MeasureList, idx) = idx > length(list)

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

@inline function snapshot{T}(list::MeasureList{T}, y::Vector{Complex{T}},
                             idx::Int, t::T)
    list.tidx_max < idx && measure_list_update(list, idx)
    snapshot(list.cur_measure, y, idx - list.tidx_offset, t - list.t_offset)
end

immutable FullMeasure{T} <: AbstractMeasure{T}
    ys::Matrix{T}
    FullMeasure(nsteps) = new(Matrix{T}(2, nsteps + 1))
    # Standard builder interface
    FullMeasure(idxs, dt) = FullMeasure{T}(length(idxs) - 1)
end

@inline function snapshot(measure::FullMeasure, y, idx, t)
    @inbounds measure.ys[1, idx] = abs2(y[1])
    @inbounds measure.ys[2, idx] = abs2(y[2])
    nothing
end

immutable SingleMeasure{T} <: AbstractMeasure{T}
    y::Vector{Complex{T}}
    # Standard builder interface
    function SingleMeasure(idxs, dt)
        @assert length(idxs) == 1
        new(Vector{Complex{T}}(2))
    end
end

@inline function snapshot(measure::SingleMeasure, y, idx, t)
    @inbounds measure.y[:] = y
    nothing
end

end
