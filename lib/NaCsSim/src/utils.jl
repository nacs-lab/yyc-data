#!/usr/bin/julia

module Utils

# Do not make this an AbstractArray since I'm not sure how useful it would be.
# In fact, this is so special purpose that I can't think of many generic
# operations that I want to have on it....
type WaveFunc{T,N} # <: AbstractArray{Complex{T},N}
    sum::T
    sz::NTuple{N,Int}
    sparse::Vector{Tuple{NTuple{N,Int},Complex{T}}}
    WaveFunc(sz::NTuple{N,Int}) =
        new(0, sz, Tuple{NTuple{N,Int},Complex{T}}[])
end
(::Type{WaveFunc{T,N}}){T,N}(sz::Vararg{Int,N}) = WaveFunc{T,N}(sz)
(::Type{WaveFunc{T}}){T,N}(sz::NTuple{N,Int}) = WaveFunc{T,N}(sz)
(::Type{WaveFunc{T}}){T,N}(sz::Vararg{Int,N}) = WaveFunc{T,N}(sz)

function zero!(ary::WaveFunc)
    ary.sum = 0
    empty!(ary.sparse)
    return
end
Base.size(ary::WaveFunc) = ary.sz

end
