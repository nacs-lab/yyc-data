#!/usr/bin/julia -f

# Some utility functions and types

module Utils

import Base: *, /, \, -, +, ==

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
    Expr(:meta, esc(x))
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

export Vec3D

"""
A 3D vector
"""
immutable Vec3D{T<:Number}
    x::T
    y::T
    z::T
end

@inline (==)(vec1::Vec3D, vec2::Vec3D) =
    vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z

@inline Base.isequal(vec1::Vec3D, vec2::Vec3D) =
    (isequal(vec1.x, vec2.x) && isequal(vec1.y, vec2.y) &&
     isequal(vec1.z, vec2.z))

"""
Sum of 3D vectors
"""
@generated function +(vec1::Vec3D, vecs::Vec3D...)
    @meta_expr inline
    len = length(vecs)
    quote
        $(Expr(:meta, :inline))
        Vec3D(+(vec1.x, $([:(vecs[$i].x) for i in 1:len]...)),
              +(vec1.y, $([:(vecs[$i].y) for i in 1:len]...)),
              +(vec1.z, $([:(vecs[$i].z) for i in 1:len]...)))
    end
end

+(vec::Vec3D) = vec

@inline -(vec::Vec3D) = Vec3D(-vec.x, -vec.y, -vec.z)

@inline -(vec1::Vec3D, vec2::Vec3D) = Vec3D(vec1.x - vec2.x, vec1.y - vec2.y,
                                            vec1.z - vec2.z)

"""
Left multiply by scalar
"""
*{T<:Number}(s::T, vec::Vec3D) = Vec3D(s * vec.x, s * vec.y, s * vec.z)

"""
Right multiply by scalar
"""
*{T<:Number}(vec::Vec3D, s::T) = Vec3D(vec.x * s, vec.y * s, vec.z * s)

"""
Divide by scalar
"""
/{T<:Number}(vec::Vec3D, s::T) = Vec3D(vec.x / s, vec.y / s, vec.z / s)

"""
Divide by scalar
"""
\{T<:Number}(s::T, vec::Vec3D) = Vec3D(s \ vec.x, s \ vec.y, s \ vec.z)

"""
Cross product of two 3D vectors
"""
Base.cross(a::Vec3D, b::Vec3D) =
    Vec3D(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)

"""
Dot product of two 3D vectors
"""
*(a::Vec3D, b::Vec3D) = a.x * b.x + a.y * b.y + a.z * b.z

@inline Base.abs2(vec::Vec3D) = abs2(vec.x) + abs2(vec.y) + abs2(vec.z)
@inline Base.abs(vec::Vec3D) = sqrt(abs2(vec))

end
