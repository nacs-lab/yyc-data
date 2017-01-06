#!/usr/bin/julia

module System

import ..Utils
import ..Samplers
import ..Setup

import ..Utils: HybridArray

# Atomic state

typealias State{T<:AbstractFloat,N} NTuple{N,HybridArray{T,3}}

function (::Type{State{T,N}}){T,N}(nx, ny, nz)
    ntuple(x->HybridArray{T,3}(nx + 1, ny + 1, nz + 1), Val{N})
end

# Initial condition

immutable ThermalInit{Idx,T<:AbstractFloat}
    nx::T
    ny::T
    nz::T
end

function (init::ThermalInit{Idx,T}){Idx,T}(atomic_state::State{T})
    for ary in atomic_state
        Utils.zero!(ary)
    end
    ary = atomic_state[Idx]
    nmax = size(ary.full)
    nx = Samplers.thermal(init.nx, nmax[1] - 1)
    ny = Samplers.thermal(init.ny, nmax[2] - 1)
    nz = Samplers.thermal(init.nz, nmax[3] - 1)
    push!(ary.sparse, ((nx + 1, ny + 1, nz + 1), Complex{T}(1)))
    return
end

# Pulses

# External state / measure

end
