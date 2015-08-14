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

call{Pol,T}(::Type{TrigCache}, trans::Transition{Pol,T}, xs) =
    TrigCache{T}(xs .* trans.k)

@generated get_transition_type{T<:Transition}(::Type{T}) =
    T.parameters[1]::TransitionType

export AtomBuilder, add_state!, add_transition!

immutable AtomBuilder{T}
    state_grps::Vector{Symbol}
    states::Vector{Tuple{Symbol,T,Int}} # Array of (name, energy, grp_id)
    transitions::Dict{NTuple{2,Int},Transition{ANY,T}} # Dict of transitions
    AtomBuilder() = new(Vector{Symbol}(),
                        Vector{Pair{Symbol,T}}(),
                        Dict{NTuple{2,Int},Transition{ANY,T}}())
end

function get_state_id(builder::AtomBuilder, _name)
    name = symbol(_name)
    get_state_id(builder, name)
end

function get_state_id(builder::AtomBuilder, name::Symbol)
    @inbounds for i in 1:length(builder.states)
        builder.states[i][1] == name && return (name, i)
    end
    name, 0
end

function get_state_grp_id(builder::AtomBuilder, _grp)
    grp = symbol(_grp)
    get_state_grp_id(builder, grp)
end

function get_state_grp_id(builder::AtomBuilder, grp::Symbol)
    @inbounds for i in 1:length(builder.state_grps)
        builder.state_grps[i] == grp && return (grp, i)
    end
    grp, 0
end

function add_state!{T}(builder::AtomBuilder{T}, _name, _grp, _energy)
    name, id = get_state_id(builder, _name)
    id == 0 || throw(ArgumentError("name $name already exist at index $id"))
    grp, gid = get_state_grp_id(builder, _grp)
    if gid == 0
        push!(builder.state_grps, grp)
        gid = length(builder.state_grps)
    end
    energy = T(_energy)
    push!(builder.states, (name, energy, gid))
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

immutable InternStates{Names,N,T,Trans,TransLevels,Grps,StateGids}
    # Names::NTuple{N,Symbol}:
    #     Names of the states
    # Trans::Type{NTuple{M,Transition{Pol::TransitionType,T}}}
    #     Types of the transitions
    # TransLevels::NTuple{M,NTuple{2,Int}}
    #     The level pairs corresponding to each transition
    # Grps::NTuple{GN,Symbol}:
    #     Names of the groups
    # StateGids::NTuple{N,Int}:
    #     Group ID of each states
    energies::NTuple{N,T}
    transitions::Trans
end

@generated num_states{T<:InternStates}(::Type{T}) = T.parameters[2]::Int
@inline num_states{T<:InternStates}(::T) = num_states(T)

@generated get_transition_types{T<:InternStates}(::Type{T}) =
    (T.parameters[4].parameters...)

@generated get_transition_pairs{T<:InternStates}(::Type{T}) = T.parameters[5]

function get_state_id(intern::InternStates, _name)
    name = symbol(_name)
    get_state_id(intern, name)
end

function get_state_id{Names}(::InternStates{Names}, name::Symbol)
    @inbounds for i in 1:length(Names)
        Names[i] == name && return (name, i)
    end
    name, 0
end

@generated get_state_gids{T<:InternStates}(::Type{T}) = T.parameters[7]

@generated get_state_names{T<:InternStates}(::Type{T}) = T.parameters[1]

function call{T}(::Type{InternStates}, builder::AtomBuilder{T})
    Names = (Symbol[state[1] for state in builder.states]...)
    energies = (T[state[2] for state in builder.states]...)
    Nstates = length(Names)
    transition_pairs = collect(builder.transitions)
    transitions = (Transition{ANY,T}[trans_pair.second
                                     for trans_pair in transition_pairs]...)
    Trans = typeof(transitions)
    TransLevels = (NTuple{2,Int}[trans_pair.first
                                 for trans_pair in transition_pairs]...)

    Grps = (Symbol[grp for grp in builder.state_grps]...)
    StateGids = (Int[state[3] for state in builder.states]...)

    InternStates{Names,Nstates,T,Trans,TransLevels,Grps,StateGids}(energies,
                                                                   transitions)
end

end
