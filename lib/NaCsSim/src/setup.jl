#!/usr/bin/julia

module Setup

import FunctionWrappers: FunctionWrapper

immutable Sequence{AS,ES,IA,I,M}
    init_atom::IA
    init::I
    pulses::Vector{FunctionWrapper{Bool,Tuple{AS,ES}}}
    measure::M
end

function run{AS,ES}(seq::Sequence{AS,ES}, atomic_state::AS, extern_state::ES)
    seq.init_atom(atomic_state)
    seq.init(atomic_state, extern_state)
    @inbounds for p in seq.pulses
        if !p(atomic_state, extern_state)
            break
        end
    end
    return seq.measure(atomic_state, extern_state)
end

immutable SeqBuilder{AS,ES,IA,I,M}
    seq::Sequence{AS,ES,IA,I,M}
    pulsemap::Dict{Any,Any}
end
function (::Type{SeqBuilder{AS,ES}}){AS,ES,IA,I,M}(init_atom::IA, init::I,
                                                   measure::M)
    seq = Sequence{AS,ES,IA,I,M}(init_atom, init,
                                 FunctionWrapper{Bool,Tuple{AS,ES}}[],
                                 measure)
    return SeqBuilder(seq, Dict{Any,Any}())
end

compile_pulse(pulse) = pulse

function add_pulse(builder::SeqBuilder, pulse)
    real_pulse = get!(builder.pulsemap, pulse) do
        compile_pulse(pulse)
    end
    pulses = builder.seq.pulses
    push!(pulses, eltype(pulses)(real_pulse))
    return
end

end
