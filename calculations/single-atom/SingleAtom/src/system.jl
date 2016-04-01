#!/usr/bin/julia -f

###
# Description of the whole system

module System

using Compat
using ..Atomic
using ..Optical
using ..Utils

import ..Atomic: add_state!, add_transition!, get_state_id, get_state_grp_id
import ..Atomic: num_states, get_transition_types, get_transition_pairs
import ..Atomic: get_state_gids, get_state_names

export AbstractPotential, HarmonicPotential, ZeroPotential
export get_potential, get_kinetic

abstract AbstractPotential{T}

"""
get_potential{T}(::AbstractPotential{T}, m::T, x::T)

Get the potential energy
"""
function get_potential end

"""
Get the kinetic energy
"""
@inline get_kinetic{T}(m::T, k::T) = k^2 / (2 * m)

immutable HarmonicPotential{T} <: AbstractPotential{T}
    ω::T
end

@inline function get_potential{T}(h::HarmonicPotential{T}, m::T, x::T)
    ω = h.ω
    sign(ω) * m * h.ω^2 * x^2 / 2
end

immutable ZeroPotential{T} <: AbstractPotential{T}
end

@inline get_potential{T}(h::ZeroPotential{T}, m::T, x::T) = T(0)

export SystemBuilder, add_state!, add_transition!, add_potential!, add_drive!

# System builder
#     An atom with potentials and drives
immutable SystemBuilder{T}
    atom::AtomBuilder{T}
    potentials::Dict{Int,AbstractPotential{T}} # Dict of potentials
    default_potential::Ref{AbstractPotential{T}}
    drives::Vector{Tuple{Union{Drive{ANY,Complex{T}},Drive{ANY,T}},Int,Int}}
    @inline SystemBuilder() = new(AtomBuilder{T}(),
                                  Dict{Int,AbstractPotential{T}}(),
                                  Ref{AbstractPotential{T}}(ZeroPotential{T}()),
                                  Vector{Drive{ANY,Union{T,Complex{T}}}}())
end

@inline add_state!(builder::SystemBuilder, args...) =
    add_state!(builder.atom, args...)

@inline add_transition!(builder::SystemBuilder, args...) =
    add_transition!(builder.atom, args...)

@inline get_state_id(builder::SystemBuilder, args...) =
    get_state_id(builder.atom, args...)

@inline get_state_grp_id(builder::SystemBuilder, args...) =
    get_state_grp_id(builder.atom, args...)

@inline function add_potential!{T}(builder::SystemBuilder{T},
                                   p::AbstractPotential{T})
    builder.default_potential[] = p
end

@inline function add_potential!{T}(builder::SystemBuilder{T},
                                   p::AbstractPotential{T}, _name)
    name, id = get_state_id(builder, _name)
    id == 0 && throw(ArgumentError("Invalid state: $_name"))
    if id in keys(builder.potentials)
        throw(ArgumentError("Potential of state $_name already set"))
    end
    builder.potentials[id] = p
    builder
end

@inline function add_drive!(builder::SystemBuilder, drive, _from_grp, _to_grp)
    from_grp, from_id = get_state_grp_id(builder, _from_grp)
    to_grp, to_id = get_state_grp_id(builder, _to_grp)
    from_id == to_id && throw(ArgumentError(string("The drive should be ",
                                                   "between different ",
                                                   "state groups")))
    from_id == 0 && throw(ArgumentError("Invalid 'from' state group: $_from_grp"))
    to_id == 0 && throw(ArgumentError("Invalid 'to' state group: $_to_grp"))
    push!(builder.drives, (drive, from_id, to_id))
    builder
end

export MotionSystem

# MotionSystem

# The atom is moving along the x-axis
immutable MotionSystem{Ax,T,PotIdxs,Intern<:InternStates,Pots,Dris,DriGrps}
    # Ax:
    #     The quantization axis of the atom
    # PotIdxs:
    #     Potential indexes into the potentials Tuple for each internal states
    # DriGrps:
    #     State groups id for each drive
    mass::T
    intern::Intern # Internal states
    potentials::Pots # Collections of potentials
    drives::Dris
end

(::Type{MotionSystem}){T}(ax::Vec3D{T}, mass, builder::SystemBuilder{T}) =
    MotionSystem{ax}(mass, builder)

@generated get_value_type{T<:MotionSystem}(::Type{T}) = T.parameters[2]

@generated get_drive_types{T<:MotionSystem}(::Type{T}) =
    (T.parameters[6].parameters...)

@generated get_potential_types{T<:MotionSystem}(::Type{T}) =
    (T.parameters[5].parameters...)

@generated num_states{T<:MotionSystem}(::Type{T}) = num_states(T.parameters[4])

@generated get_transition_types{T<:MotionSystem}(::Type{T}) =
    get_transition_types(T.parameters[4])

@generated get_state_gids{T<:MotionSystem}(::Type{T}) =
    get_state_gids(T.parameters[4])

@generated get_state_names{T<:MotionSystem}(::Type{T}) =
    get_state_names(T.parameters[4])

@generated get_transition_pairs{T<:MotionSystem}(::Type{T}) =
    get_transition_pairs(T.parameters[4])

@generated get_potential_idxs{T<:MotionSystem}(::Type{T}) = T.parameters[3]

@generated get_quant_axis{T<:MotionSystem}(::Type{T}) = T.parameters[1]

@inline get_state_id(sys::MotionSystem, args...) =
    get_state_id(sys.intern, args...)

@generated get_drive_gids{T<:MotionSystem}(::Type{T}) = T.parameters[7]


function (::Type{MotionSystem{Ax}}){Ax,T}(_mass, builder::SystemBuilder{T})
    Base.typeassert(Ax, Vec3D{T})
    mass = T(_mass)
    intern = InternStates(builder.atom)
    Intern = typeof(intern)
    N = num_states(intern)

    _pots = Vector{AbstractPotential{T}}()
    _pot_map = Dict{AbstractPotential{T},Int}()
    _pot_idxs = Vector{Int}(N)
    for i in 1:N
        pot = if i in keys(builder.potentials)
            builder.potentials[i]
        else
            builder.default_potential[]
        end

        if !(pot in keys(_pot_map))
            push!(_pots, pot)
            _pot_map[pot] = length(_pots)
        end
        _pot_idxs[i] = _pot_map[pot]
    end
    potentials = (_pots...)
    Pots = typeof(potentials)
    drives = ([drive[1] for drive in builder.drives]...)
    Dris = typeof(drives)
    PotIdxs = (_pot_idxs...)
    DriGrps = (NTuple{2,Int}[(drive[2], drive[3])
                             for drive in builder.drives]...)
    MotionSystem{Ax,T,PotIdxs,Intern,Pots,Dris,DriGrps}(mass, intern,
                                                        potentials, drives)
end

end
