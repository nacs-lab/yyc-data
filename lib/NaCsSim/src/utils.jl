#!/usr/bin/julia

module Utils

# Do not make this an AbstractArray since I'm not sure how useful it would be.
# In fact, this is so special purpose that I can't think of many generic
# operations that I want to have on it....
type HybridArray{T,N} # <: AbstractArray{Complex{T},N}
    isfull::Bool
    sum::T
    full::Array{Complex{T},N}
    sparse::Vector{Tuple{NTuple{N,Int},Complex{T}}}
    HybridArray(sz::NTuple{N,Int}) =
        new(false, 0, Array{Complex{T},N}(sz), Tuple{NTuple{N,Int},Complex{T}}[])
end
(::Type{HybridArray{T,N}}){T,N}(sz::Vararg{Int,N}) = HybridArray{T,N}(sz)
(::Type{HybridArray{T}}){T,N}(sz::NTuple{N,Int}) = HybridArray{T,N}(sz)
(::Type{HybridArray{T}}){T,N}(sz::Vararg{Int,N}) = HybridArray{T,N}(sz)

function zero!(ary::HybridArray)
    ary.isfull = false
    ary.sum = 0
    empty!(ary.sparse)
    return
end

end
