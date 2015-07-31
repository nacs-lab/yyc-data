#!/usr/bin/julia -f

# Some utility functions and types

module Utils

# Structure of Arrays

# Helper types and functions to help using them with complex numbers
# and with array of structures format.

using StructsOfArrays

export StructOfArrays, SoCArray, SoCVector, SoCMatrix

typealias SoCArray{T,N} StructOfArrays{Complex{T},N,NTuple{2,Array{T,N}}}
typealias SoCVector{T} SoCArray{T,1}
typealias SoCMatrix{T} SoCArray{T,2}

@inline function Base.unsafe_copy!{T,N}(dest::SoCArray{T,N},
                                        src::Array{Complex{T},N})
    len = length(src)
    BLAS.blascopy!(len, Ptr{T}(pointer(src)), 2,
                   pointer(dest.arrays[1]), 1)
    BLAS.blascopy!(len, Ptr{T}(pointer(src)) + sizeof(T), 2,
                   pointer(dest.arrays[2]), 1)
    dest
end

@inline function Base.unsafe_copy!{T,N}(dest::Array{Complex{T},N},
                                        src::SoCArray{T,N})
    len = length(src)
    BLAS.blascopy!(len, pointer(src.arrays[1]), 1,
                   Ptr{T}(pointer(dest)), 2)
    BLAS.blascopy!(len, pointer(src.arrays[2]), 1,
                   Ptr{T}(pointer(dest)) + sizeof(T), 2)
    dest
end

export TrigCache

"""
Store pre-compute values of sin's and cos's on the grid points
"""
immutable TrigCache{T}
    sins::Vector{T}
    coss::Vector{T}
    function TrigCache(θs)
        nele = length(θs)
        sins = Vector{T}(nele)
        coss = Vector{T}(nele)
        @inbounds for i in 1:nele
            sins[i] = sin(T(θs[i]))
            coss[i] = cos(T(θs[i]))
        end
        new(sins, coss)
    end
end

TrigCache{T}(θs::AbstractArray{T}) = TrigCache{T}(θs)

export @meta_expr

macro meta_expr(x)
    Expr(:meta, x)
end

export sum2average

"""
Convert from Σ(x) and Σ(x²) to the average and uncertainty of average
"""
function sum2average(s, s2, count)
    avg = s / count
    avg2 = s2 / count
    std = (avg2 - avg^2) / (count - 1)
    # rounding errors can make small std smaller than zero
    unc = std <= 0 ? zero(std) : sqrt(std)
    avg, unc
end

"""
Cross product of two vectors
"""
Base.cross{T1<:Number,T2<:Number}(a::NTuple{3,T1}, b::NTuple{3,T2}) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
     a[1] * b[2] - a[2] * b[1])

end
