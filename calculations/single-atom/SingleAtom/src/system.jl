#!/usr/bin/julia -f

###
# Description of the whole system

module System

using ..Atomic
using ..Optical
using ..Utils

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
@inline function get_kinetic{T}(m::T, k::T)
    k^2 / (2 * m)
end

immutable HarmonicPotential{T} <: AbstractPotential{T}
    ω::T
end

@inline get_potential{T}(h::HarmonicPotential{T}, m::T, x::T) =
    m * h.ω^2 * x^2 / 2

immutable ZeroPotential{T} <: AbstractPotential{T}
end

@inline get_potential{T}(h::HarmonicPotential{T}, m::T, x::T) = T(0)

export SystemBuilder, add_state!, add_transition!, add_potential!, add_drive!

# System builder
#     An atom with potentials and drives
immutable SystemBuilder{T}
    atom::AtomBuilder{T}
    potentials::Dict{Int,AbstractPotential{T}} # Dict of potentials
    default_potential::Ref{AbstractPotential{T}}
    drives::Vector{Union{Drive{ANY,Complex{T}},Drive{ANY,T}}}
    @inline SystemBuilder() = new(AtomBuilder{T}(),
                                  Dict{Int,AbstractPotential{T}}(),
                                  Ref{AbstractPotential{T}}(ZeroPotential{T}()),
                                  Vector{Drive{ANY,Union{T,Complex{T}}}}())
end

@inline Atomic.add_state!(builder::SystemBuilder, args...) =
    Atomic.add_state!(builder.atom, args...)

@inline Atomic.add_transition!(builder::SystemBuilder, args...) =
    Atomic.add_transition!(builder.atom, args...)

@inline Atomic.get_state_id(builder::SystemBuilder, args...) =
    Atomic.get_state_id(builder.atom, args...)

@inline function add_potential!{T}(builder::SystemBuilder{T},
                                   p::AbstractPotential{T})
    builder.default_potential[] = p
end

@inline function add_potential!{T}(builder::SystemBuilder{T},
                                   p::AbstractPotential{T}, _name)
    name, id = Atomic.get_state_id(builder, _name)
    id == 0 && throw(ArgumentError("Invalid state: $_name"))
    if id in keys(builder.potentials)
        throw(ArgumentError("Potential of state $_name already set"))
    end
    builder.potentials[id] = p
    builder
end

@inline function add_drive!(builder::SystemBuilder, drive)
    push!(builder.drives, drive)
    builder
end

export MotionSystem

# MotionSystem

immutable MotionSystem{Ax,T,PotIdxs,Intern<:InternStates,Pots,Dris}
    mass::T
    intern::Intern
    potentials::Pots
    drives::Dris
end

call{T}(::Type{MotionSystem}, ax::Vec3D{T}, mass, builder::SystemBuilder{T}) =
    MotionSystem{ax}(mass, builder)

function call{Ax,T}(::Type{MotionSystem{Ax}}, _mass, builder::SystemBuilder{T})
    Base.typeassert(Ax, Vec3D{T})
    mass = T(_mass)
    intern = InternStates(builder.atom)
    Intern = typeof(intern)
    N = Atomic.num_states(intern)

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
    drives = (builder.drives...)
    Dris = typeof(drives)
    PotIdxs = (_pot_idxs...)
    MotionSystem{Ax,T,PotIdxs,Intern,Pots,Dris}(mass, intern,
                                                potentials, drives)
end

end
