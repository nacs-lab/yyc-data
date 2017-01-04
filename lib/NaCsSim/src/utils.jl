#!/usr/bin/julia

module Utils

using Base.Cartesian
import NaCsCalc: Trap

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

function sample_sideband(n::Int, η, nmax::Int)
    # Estimate the range of final states with none zero matrix element.
    # By starting with the ones with high probability we can minimize the
    # expectation value of evaluation time.

    # I can't find any good bound on the matrix element.
    # The following is based on a combination of observation and guessing.
    # Since we loop the whole range in the end, it doesn't affect correctness
    # if the estimation is off.

    # We estimate the range based on the max k values.
    # We effectively assumes that most population of a states are at the max k
    # i.e. `±k_max`. These are shifted to `k±k_max` after a kick `k`. The range
    # of states with high overlap is then determined by overlapping the max k
    # of the final states to these two values, i.e. `k′_max = |k±k_max|`
    # Max k for state n
    # √(2n + 1) * √(mω / ħ)
    # Displacement in k
    # √(2) * η * √(mω / ħ)
    # k / √(mω / ħ) range
    # |√(2n + 1) - √(2) * η| ~ (√(2n + 1) + √(2) * η)
    # State with significant overlap
    # n + η² ± √(4n + 2)η

    n_center = n + η^2
    n_width = √(4n + 2) * η
    n_start = max(0, floor(Int, n_center - n_width))
    n_end = min(ceil(Int, n_center + n_width), nmax)

    v::Float64 = rand()

    for i in n_start:n_end
        v -= Trap.sideband(n, i, η)^2
        v <= 0 && return i
    end
    nleft1 = n_start
    nleft2 = nmax - n_end
    for i in 1:min(nleft1, nleft2)
        i1 = n_start - i
        v -= Trap.sideband(n, i1, η)^2
        v <= 0 && return i1
        i2 = n_end + i
        v -= Trap.sideband(n, i2, η)^2
        v <= 0 && return i2
    end
    if nleft1 > nleft2
        for i in (n_start - nleft2 - 1):-1:0
            v -= Trap.sideband(n, i, η)^2
            v <= 0 && return i
        end
    elseif nleft2 > nleft1
        for i in (n_end + nleft1 + 1):nmax
            v -= Trap.sideband(n, i, η)^2
            v <= 0 && return i
        end
    end
    return -1
end

end
