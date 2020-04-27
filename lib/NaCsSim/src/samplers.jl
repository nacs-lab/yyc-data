#!/usr/bin/julia

module Samplers

using Base.Cartesian
import NaCsCalc: Trap

@generated default_index(::Val{N}, ::Val{M}=Val{0}()) where {N,M} = ntuple(i->M, N)

function sideband(n::Int, η::T, nmax::Int, rng) where {T}
    # Estimate the range of final states with non-zero matrix elements.
    # By starting with the ones with high probability we can minimize the
    # average evaluation time.

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

    v::T = rand(rng)

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

function emission(::Type{T}, isσ::Bool, rng) where T<:AbstractFloat
    # Returns `(cosθ, φ)`. The caller should be ready to handle `|cosθ| > 1`
    # The PDFs of the `θ` distribution are
    # `3 / 4 * (1 - cos²θ) * sinθ` for π light and
    # `3 / 8 * (1 + cos²θ) * sinθ` for σ± light
    # The corresponding CDFs are
    # `1 / 4 * (3cosθ - cos³θ) + 1 / 2` for π light and
    # `1 / 8 * (3cosθ + cos³θ) + 1 / 2` for σ± light
    # For a random number `v` we need to solve
    # `3cosθ - cos³θ = 4v - 2` for π light and
    # `3cosθ + cos³θ = 8v - 4` for σ± light
    # The (real) solution is
    # `cosθ = x + 1 / x` where `x = ∛(2√(v² - v) - 2v + 1)` for π light and
    # `cosθ = x - 1 / x` where `x = ∛(√(16v² - 16v + 5) + 4v - 2)` for σ± light
    # Note that the `x` for π polarization is complex
    v = T(rand(rng)) # This is faster than `rand(T)`....
    φ = T(rand(rng)) * T(2π)
    if isσ
        y = muladd(v, 4, -2)
        x = @fastmath cbrt(sqrt(muladd(y, y, 1)) + y)
        return x - 1 / x, φ
    else
        θ′ = @fastmath acos(muladd(T(-2), v, T(1)))
        cosθ = @fastmath cos(muladd(θ′, T(1 / 3), - T(2π / 3))) * 2
        return cosθ, φ
    end
end

unit_3d_transformation(::Type{T}) where T = ((T(1), T(0), T(0)),
                                             (T(0), T(1), T(0)),
                                             (T(0), T(0), T(1)))

vec3d_norm2(vec::NTuple{3,<:Number}) = abs2(vec[1]) + abs2(vec[2]) + abs2(vec[3])
vec3d_cross(vec1::NTuple{3,<:Number}, vec2::NTuple{3,<:Number}) =
    (muladd(vec1[2], vec2[3], -(vec1[3] * vec2[2])),
     muladd(vec1[3], vec2[1], -(vec1[1] * vec2[3])),
     muladd(vec1[1], vec2[2], -(vec1[2] * vec2[1])))

"""
    vec3d_to_trans(ax::NTuple{3,T}) where T<:AbstractFloat

Given an (unnormlized) 3D vector `ax` compute 3 orthonormal vectors in which the first
vector is parallel with `ax`. The direction of the two other vectors are arbitrarily chosen.
"""
function vec3d_to_trans(ax::NTuple{3,T}) where T<:AbstractFloat
    na = vec3d_norm2(ax)
    @assert na != 0
    ax = ax ./ sqrt(na)
    v1 = vec3d_cross((1, 0, 0), ax)
    v2 = vec3d_cross((0, 1, 0), ax)
    v3 = vec3d_cross((0, 0, 1), ax)
    n1 = vec3d_norm2(v1)
    n2 = vec3d_norm2(v2)
    n3 = vec3d_norm2(v3)
    if n3 >= n1
        if n3 >= n2
            v = v3
            nv = n3
        else
            # n3 >= n1
            # n2 > n3
            v = v2
            nv = n2
        end
    elseif n1 >= n2
        v = v1
        nv = n1
    else
        # n1 > n3
        # n2 > n1
        v = v2
        nv = n2
    end
    v = v ./ sqrt(nv)
    # Should be normalized already
    vy = vec3d_cross(ax, v)
    return (ax, v, vy)
end

# This currently can't handle misalignment between the quantization axis
# and trap axis. Hopefully it's not very important
function op(n_init::NTuple{3,Int}, n_max::NTuple{3,Int}, ηs::NTuple{3,T}, ηdri::NTuple{3,T},
            isσ::Bool, rng,
            trans::NTuple{3,NTuple{3,T}}=unit_3d_transformation(T)) where T<:AbstractFloat
    cosθ, φ = emission(T, isσ, rng)
    η1 = cosθ
    if -1 < cosθ < 1
        # @fastmath on comparison is currently problematic
        sinθ = @fastmath sqrt(1 - cosθ * cosθ)
        s, c = @fastmath sincos(φ)
        η2 = sinθ * c
        η3 = sinθ * s
    else
        η2 = zero(T)
        η3 = zero(T)
    end
    transx, transy, transz = trans
    ηx = ηs[1] * muladd(η1, transx[1], muladd(η2, transy[1], η3 * transz[1]))
    ηy = ηs[2] * muladd(η1, transx[2], muladd(η2, transy[2], η3 * transz[2]))
    ηz = ηs[3] * muladd(η1, transx[3], muladd(η2, transy[3], η3 * transz[3]))
    ηx = abs(ηx - ηdri[1])
    ηy = abs(ηy - ηdri[2])
    ηz = abs(ηz - ηdri[3])
    nx = sideband(n_init[1], ηx, n_max[1], rng)
    nx == -1 && @goto escape
    ny = sideband(n_init[2], ηy, n_max[2], rng)
    ny == -1 && @goto escape
    nz = sideband(n_init[3], ηz, n_max[3], rng)
    nz == -1 && @goto escape
    return (nx, ny, nz)
    @label escape
    return default_index(Val{3}(), Val{-1}())
end

function thermal(nbar, nmax, rng)
    α = @fastmath 1 / log(nbar / (nbar + 1))
    while true
        n = @fastmath floor(Int, log(rand(rng)) * α)
        if n <= nmax
            return n
        end
    end
end

function decay(rates::AbstractArray{T}, weights::AbstractArray{T}, rng) where {T}
    total = sum(weights)
    v = T(rand(rng)) * total
    nstates = length(rates)
    if length(weights) != nstates
        throw(ArgumentError("Decay rates and weights should have the same size"))
    end
    @inbounds for i in 1:nstates
        w = weights[i]
        old_v = v
        v -= w
        @fastmath if v < 0
            return -log(old_v / w) / rates[i], i
        end
    end
    return T(Inf), 1
end

decay(rate::T, rng) where {T<:AbstractFloat} = -log(T(rand(rng))) / rate
function select(total::T, weights::Union{AbstractArray{T},Tuple{Vararg{T}}}, rng) where T
    v = T(rand(rng)) * total
    @inbounds for i in 1:length(weights)
        v -= weights[i]
        v <= 0 && return i
    end
    return 1
end

end
