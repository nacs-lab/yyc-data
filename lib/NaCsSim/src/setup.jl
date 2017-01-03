#!/usr/bin/julia

module Setup

import FunctionWrappers: FunctionWrapper

immutable Sequence{AS,ES,I,M}
    init::I
    pulses::Vector{FunctionWrapper{Void,Tuple{AS,ES}}}
    measure::M
end

function run{AS,ES}(seq::Sequence, atomic_state::AS, extern_state::ES)
    seq.init(atomic_state, extern_state)
    @inbounds for p in seq.pulses
        p(atomic_state, extern_state)
    end
    return seq.measure(atomic_state, extern_state)
end

end
