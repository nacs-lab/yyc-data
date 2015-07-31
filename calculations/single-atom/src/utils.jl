#!/usr/bin/julia -f

# Some utility functions and types

module Utils

using StructsOfArrays

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

macro meta_expr(x)
    Expr(:meta, x)
end

function sum2average(s, s2, count)
    avg = s / count
    avg2 = s2 / count
    std = (avg2 - avg^2) / (count - 1)
    # rounding errors can make small std smaller than zero
    unc = std <= 0 ? zero(std) : sqrt(std)
    avg, unc
end

end
