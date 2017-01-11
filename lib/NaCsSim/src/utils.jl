#!/usr/bin/julia

module Utils

immutable WaveFunc{T,N}
    abs2::Base.RefValue{T}
    sz::NTuple{N,Int}
    sparse::Vector{Tuple{NTuple{N,Int},Complex{T}}}
    WaveFunc(sz::NTuple{N,Int}) =
        new(Ref{T}(0), sz, Tuple{NTuple{N,Int},Complex{T}}[])
end
(::Type{WaveFunc{T,N}}){T,N}(sz::Vararg{Int,N}) = WaveFunc{T,N}(sz)
(::Type{WaveFunc{T}}){T,N}(sz::NTuple{N,Int}) = WaveFunc{T,N}(sz)
(::Type{WaveFunc{T}}){T,N}(sz::Vararg{Int,N}) = WaveFunc{T,N}(sz)

function zero!(ary::WaveFunc)
    ary.abs2[] = 0
    empty!(ary.sparse)
    return
end
Base.size(ary::WaveFunc) = ary.sz
@inline Base.abs2{T}(ary::WaveFunc{T}) = ary.abs2[]
function Base.push!{T}(ary::WaveFunc{T}, ns, v)
    v2 = abs2(v)
    ary.abs2[] += v2
    push!(ary.sparse, (ns, Complex{T}(v)))
end

end
