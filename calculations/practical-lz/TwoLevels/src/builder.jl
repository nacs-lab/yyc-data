#!/usr/bin/julia -f

module Builders

import ..TwoLevels: AbstractDrive, AbstractMeasure, MeasureList, Sequence

immutable DriveWrapper{T}
    len::Int
    drive::AbstractDrive{T}
end

type MeasureWrapper
    i_start::Int
    i_end::Int
    measure
end
# MeasureList expect ((imin, imax), m)
Base.start(m::MeasureWrapper) = 1
Base.next(m::MeasureWrapper, idx) = if idx == 1
    (m.i_start, m.i_end), 2
elseif idx == 2
    m.measure, 3
else
    throw(BoundsError(m, idx))
end
Base.done(m::MeasureWrapper, idx) = idx > 2

type SeqBuilder{T}
    dt::T # Current time, starts with T(0)
    cur_t::Int # Current point index, starts with 1
    drives::Vector{DriveWrapper{T}}
    measures::Vector{MeasureWrapper}
    cur_measure_t::Int # starts with 1 (t = 0), initialized to 0
    cur_measure
    SeqBuilder(dt) = new(dt, 1, DriveWrapper{T}[], MeasureWrapper[], 0, nothing)
end
SeqBuilder{T}(dt::T) = SeqBuilder{T}(dt)

function add_drive{T}(builder::SeqBuilder{T}, t, drive::AbstractDrive)
    t_i = round(Int, t / builder.dt)
    t_i >= 1 || throw(ArgumentError("Drive length ($t) is too short"))
    builder.cur_t += t_i
    push!(builder.drives, DriveWrapper{T}(t_i, drive))
    nothing
end

function start_measure{T}(builder::SeqBuilder{T}, t::Number, measure...)
    if builder.cur_measure !== nothing
        throw(ArgumentError("There is already a on-going measurement"))
    end
    # plus one so that zero time delay is the next point
    t_i = round(Int, t / builder.dt) + builder.cur_t
    if t_i == builder.cur_measure_t
        t_i += 1
    elseif t_i < builder.cur_measure_t
        throw(ArgumentError("Overlap with previous measurement"))
    end
    builder.cur_measure_t = t_i
    if length(measure) == 1 && isa(measure[1], AbstractMeasure)
        measure = measure[1]
    end
    builder.cur_measure = measure
    nothing
end

start_measure{T}(builder::SeqBuilder{T}, measure...) =
    start_measure(builder, zero(T), measure...)

function finish_measure{T}(builder::SeqBuilder{T}, t::Number=zero(T))
    if builder.cur_measure === nothing
        throw(ArgumentError("There is no on-going measurement to finish"))
    end
    t_i = round(Int, t / builder.dt) + builder.cur_t
    if t_i < builder.cur_measure_t
        throw(ArgumentError("Negative measurement length"))
    end
    push!(builder.measures, MeasureWrapper(builder.cur_measure_t, t_i,
                                           builder.cur_measure))
    builder.cur_measure_t = t_i
    builder.cur_measure = nothing
    nothing
end

function build{T}(builder::SeqBuilder{T})
    t_len = builder.cur_t
    if builder.cur_measure !== nothing
        throw(ArgumentError("Measure need to be explicitely finished"))
    elseif builder.cur_measure_t > t_len
        throw(ArgumentError("Measure longer than sequence"))
    end
    dt = builder.dt
    measure = MeasureList{T}(1:t_len, dt, builder.measures)
    drives = Pair{Int}[(d.len=>d.drive) for d in builder.drives]
    Sequence{T}(dt, t_len, (drives...), measure)
end

end
