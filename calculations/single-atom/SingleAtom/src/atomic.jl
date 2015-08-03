#!/usr/bin/julia -f

module Atomic

using ..Utils
import ..Optical
import Base: *

export TransitionType, Trans_σ⁺, Trans_σ⁻, Trans_π

@enum TransitionType Trans_σ⁺ Trans_σ⁻ Trans_π

"""
The type of an optical dipole transition. σ⁺, σ⁻, or π
end JuliaLang/julia#12419
"""
TransitionType

"""
Overlap of a drive with a dipole transition
"""
function *(ax_trans::Tuple{Vec3D,TransitionType}, amp::Vec3D)
    # A = a * π + b * σ⁺ + c * σ⁻
    ax = ax_trans[1]
    trans = ax_trans[2]
    ax_scale = 1 / abs(ax)
    norm_ax = ax * ax_scale
    if trans == Trans_π
        # ax ⋅ π = 1
        # ax ⋅ σ⁺ = 0
        # ax ⋅ σ⁻ = 0
        return abs(norm_ax * amp)
    end
    # TODO figure out the sign of σ⁺ and σ⁻
    # ax × π = 0
    # ax × σ⁺ = i σ⁺
    # ax × σ⁻ = -i σ⁻
    trans_diff = norm_ax × amp # i (σ⁺ - σ⁻)
    # ax × (ax × π) = 0
    # ax × (ax × σ⁺) = -σ⁺
    # ax × (ax × σ⁻) = -σ⁻
    trans_amp = norm_ax × trans_diff # -(σ⁺ + σ⁻)
    σ_sign = ifelse(trans == Trans_σ⁺, -1im, 1im)
    return abs(trans_diff + σ_sign * trans_amp) / 2
end

*(amp::Vec3D, ax_trans::Tuple{Vec3D,TransitionType}) = ax_trans * amp

export Transition

"""
An optical dipole transition
"""
immutable Transition{Pol,T} # Pol::TransitionType
    k::T # wave vector (projected on the axis)
    Γ::T # Line width / decay rate
    α::T # Dipole moment
end

@inline call{Pol,T}(::Type{Transition{Pol}}, args::T...) =
    Transition{Pol,T}(args...)

export AtomBuilder, add_state!, add_transition!

immutable AtomBuilder{T}
    states::Vector{Pair{Symbol,T}} # Array of name=>energy
    transitions::Dict{NTuple{2,Int},Transition{ANY,T}} # Dict of transitions
    AtomBuilder() = new(Vector{Pair{Symbol,T}}(),
                        Dict{NTuple{2,Int},Transition{ANY,T}}())
end

function get_state_id(builder::AtomBuilder, _name)
    name = symbol(_name)
    name, get_state_id(builder, name)
end

function get_state_id(builder::AtomBuilder, name::Symbol)
    @inbounds for i in 1:length(builder.states)
        builder.states[i].first == name && return i
    end
    return 0
end

function add_state!{T}(builder::AtomBuilder{T}, _name, _energy)
    name, id = get_state_id(builder, _name)
    id == 0 || throw(ArgumentError("name $name already exist at index $id"))
    energy = T(_energy)
    push!(builder.states, name=>energy)
    builder
end

function add_transition!{Pol,T}(builder::AtomBuilder{T}, _from, _to,
                                transition::Transition{Pol,T})
    from, from_i = get_state_id(builder, _from)
    to, to_i = get_state_id(builder, _to)
    from == to && throw(ArgumentError(string("The transition should be ",
                                             "between different states")))
    from_i == 0 && throw(ArgumentError("Invalid 'from' state: $_from"))
    to_i == 0 && throw(ArgumentError("Invalid 'to' state: $_to"))
    d = builder.transitions
    if (from_i, to_i) in keys(d) || (to_i, from_i) in keys(d)
        throw(ArgumentError("Transition $_to => $_from already exist"))
    end
    d[(from_i, to_i)] = transition
    builder
end

export InternStates

immutable InternStates{Names,N,T,Trans,TransLevels}
    # Names::NTuple{N,Symbol}:
    #     Names of the states
    # Trans::Type{NTuple{M,Transition{Pol::TransitionType,T}}}
    #     Types of the transitions
    # TransLevels::NTuple{M,NTuple{2,Int}}
    #     The level pairs corresponding to each transition
    energies::NTuple{N,T}
    transitions::Trans
end

@generated num_states{T<:InternStates}(::Type{T}) = T.parameters[2]::Int
@inline num_states{T<:InternStates}(::T) = num_states(T)

@generated get_transition_types{T<:InternStates}(::Type{T}) =
    (T.parameters[4].parameters...)

function call{T}(::Type{InternStates}, builder::AtomBuilder{T})
    Names = (Symbol[state.first for state in builder.states]...)
    energies = (T[state.second for state in builder.states]...)
    Nstates = length(Names)
    transition_pairs = collect(builder.transitions)
    transitions = (Transition{ANY,T}[trans_pair.second
                                     for trans_pair in transition_pairs]...)
    Trans = typeof(transitions)
    TransLevels = (NTuple{2,Int}[trans_pair.first
                                 for trans_pair in transition_pairs]...)
    InternStates{Names,Nstates,T,Trans,TransLevels}(energies, transitions)
end

end
