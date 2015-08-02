#!/usr/bin/julia -f

module SingleAtom

include("utils.jl")
using .Utils
export StructOfArrays, SoCArray, SoCVector, SoCMatrix
export sum2average
export Vec3D

include("optical.jl")
using .Optical

include("atomic.jl")
using .Atomic
export TransitionType, Trans_σ⁺, Trans_σ⁻, Trans_π
export Transition
export AtomBuilder, add_state!, add_transition!
export InternStates

end
