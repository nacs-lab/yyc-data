#!/usr/bin/julia

module Setup

import FunctionWrappers: FunctionWrapper

struct Sequence{AS,ES,IA,I,M}
    init_atom::IA
    init::I
    pulses::Vector{FunctionWrapper{Bool,Tuple{AS,ES}}}
    measure::M
end

function run(seq::Sequence{AS,ES}, atomic_state::AS, extern_state::ES,
             res=create_measure(seq)) where {AS,ES}
    seq.init_atom(atomic_state)
    seq.init(atomic_state, extern_state)
    @inbounds for p in seq.pulses
        if !p(atomic_state, extern_state)
            break
        end
    end
    return seq.measure(res, atomic_state, extern_state)
end

# TODO, run in parallel/run multiple sequences
function run(seq::Sequence{AS,ES}, atomic_state::AS, extern_state::ES, n::Integer) where {AS,ES}
    res = create_measure(seq)
    for i in 1:n
        run(seq, atomic_state, extern_state, res)
        if abort_measure(seq, res, i)
            return finalize_measure(seq, res, i)
        end
    end
    return finalize_measure(seq, res, n)
end

struct SeqBuilder{AS,ES,IA,I,M}
    seq::Sequence{AS,ES,IA,I,M}
    pulsemap::Dict{Any,Any}
    datacache::Dict{Any,Any}
end
function (::Type{SeqBuilder{AS,ES}})(init_atom::IA, init::I,
                                     measure::M) where {AS,ES,IA,I,M}
    seq = Sequence{AS,ES,IA,I,M}(init_atom, init,
                                 FunctionWrapper{Bool,Tuple{AS,ES}}[],
                                 measure)
    return SeqBuilder(seq, Dict{Any,Any}(), Dict{Any,Any}())
end

compile_pulse(pulse, cache) = pulse

function add_pulse(builder::SeqBuilder, pulse)
    real_pulse = get!(builder.pulsemap, pulse) do
        compile_pulse(pulse, builder.datacache)
    end
    pulses = builder.seq.pulses
    push!(pulses, eltype(pulses)(real_pulse))
    return
end

struct Dummy
end
(::Dummy)(a_s, e_s) = true
(::Dummy)(res, a_s, e_s) = nothing

struct CombinedMeasure{T<:Tuple}
    measures::T
    # Kill default constructors
    CombinedMeasure{T}(measures) where T = new(measures)
end
CombinedMeasure(measures...) = CombinedMeasure{typeof(measures)}(measures)
@generated function create_measure(m::CombinedMeasure{T}, seq) where T
    N = length(T.parameters)
    quote
        ms = m.measures
        ($((:(create_measure(ms[$i], seq)) for i in 1:N)...),)
    end
end
@generated function (m::CombinedMeasure{T})(res, state, extern_state) where {T}
    N = length(T.parameters)
    quote
        ms = m.measures
        ($((:(ms[$i](res[$i], state, extern_state)) for i in 1:N)...),)
    end
end
@generated function finalize_measure(m::CombinedMeasure{T}, res, n) where {T}
    N = length(T.parameters)
    quote
        ms = m.measures
        ($((:(finalize_measure(ms[$i], res[$i], n)) for i in 1:N)...),)
    end
end
@generated function abort_measure(m::CombinedMeasure{T}, res, n) where T
    N = length(T.parameters)
    val = true
    for i in N:-1:1
        val = :(abort_measure(ms[$i], res[$i], n) && $val)
    end
    quote
        ms = m.measures
        $val
    end
end

abort_measure(measure, res, n) = false

create_measure(seq::Sequence) = create_measure(seq.measure, seq)
finalize_measure(seq::Sequence, m, n) =
    finalize_measure(seq.measure, m, n)
abort_measure(seq::Sequence, res, n) = abort_measure(seq.measure, res, n)
create_measure(::Dummy, seq) = nothing
finalize_measure(::Dummy, m, n) = nothing

end
