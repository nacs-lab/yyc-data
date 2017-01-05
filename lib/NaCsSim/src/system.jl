#!/usr/bin/julia

module System

import ..Utils
import ..Setup

import ..Utils: HybridArray

# Atomic state

typealias State{T<:AbstractFloat,N} NTuple{N,HybridArray{T,3}}

function (::Type{State{T,N}}){T,N}(nx, ny, nz)
    ntuple(x->HybridArray{T,3}(nx, ny, nz), Val{N})
end

# Initial condition

# Pulses

# External state / measure

end
