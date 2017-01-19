#!/usr/bin/julia

module Setup

import FunctionWrappers: FunctionWrapper

immutable Sequence{AS,ES,IA,I,M}
    init_atom::IA
    init::I
    pulses::Vector{FunctionWrapper{Bool,Tuple{AS,ES}}}
    measure::M
end

function run{AS,ES}(seq::Sequence{AS,ES}, atomic_state::AS, extern_state::ES,
                    res=create_measure(seq))
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
function run{AS,ES}(seq::Sequence{AS,ES},
                    atomic_state::AS, extern_state::ES, n::Integer)
    @assert n >= 1
    res = create_measure(seq)
    run(seq, atomic_state, extern_state, res)
    for i in 2:n
        run(seq, atomic_state, extern_state, res)
    end
    return finalize_measure(seq, res, n)
end

immutable SeqBuilder{AS,ES,IA,I,M}
    seq::Sequence{AS,ES,IA,I,M}
    pulsemap::Dict{Any,Any}
    datacache::Dict{Any,Any}
end
function (::Type{SeqBuilder{AS,ES}}){AS,ES,IA,I,M}(init_atom::IA, init::I,
                                                   measure::M)
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

immutable Dummy
end
(::Dummy)(a_s, e_s) = true
(::Dummy)(res, a_s, e_s) = nothing

create_measure(seq::Sequence) = create_measure(seq.measure, seq)
finalize_measure(seq::Sequence, m, n) =
    finalize_measure(seq.measure, m, n)
create_measure(::Dummy, seq) = nothing
finalize_measure(::Dummy, m, n) = nothing

end
