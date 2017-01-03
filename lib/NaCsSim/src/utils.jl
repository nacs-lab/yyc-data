#!/usr/bin/julia

module Utils

using Base.Cartesian

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

@generated default_index{N}(::Val{N}) = ntuple(i->0, N)

@generated function sample{T,N}(ary::Array{T,N}, thresh)
    quote
        $(Expr(:meta, :inline))
        @inbounds @nloops $N i ary begin
            thresh -= abs2(@nref $N ary i)
            thresh <= 0 && return @ntuple $N i
        end
        return default_index($(Val{N}()))
    end
end

function sample{T,N}(ary::HybridArray{T,N})
    thresh = rand(T) * ary.sum
    @inbounds if ary.isfull
        return sample(ary.full, thresh)
    else
        sparse_ary = ary.sparse
        for (idx, v) in sparse_ary
            thresh -= abs2(v)
            if thresh <= 0
                return idx
            end
        end
        return default_index(Val{N}())
    end
end

end
