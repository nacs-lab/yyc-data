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

end
