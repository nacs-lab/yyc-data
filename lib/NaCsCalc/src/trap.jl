#!/usr/bin/julia

module Trap

import ..Utils
using SpecialFunctions

function η(m::T1, freq::T2, k::T3) where {T1,T2,T3}
    T = float(promote_type(T1, T2, T3))
    ħ = T(1.0545718e-34)
    z_0 = √(ħ / 2 / m / (T(2π) * freq))
    z_0 * k
end

function (sideband(n1::Integer, n2::Integer, η::T)::T) where {T<:AbstractFloat}
    if n1 < 0 || n2 < 0
        return 0
    elseif η == 0
        if n1 == n2
            return 1
        else
            return 0
        end
    end
    # Ref http://journals.aps.org/pra/pdf/10.1103/PhysRevA.20.1521
    # Δn ≡ |n1 - n2|
    # n₋ ≡ min(n1, n2)
    # n₊ ≡ max(n1, n2)
    #   ⟨n1|exp(ikx)|n2⟩
    # = ⟨n1|exp(iη(a + a†))|n2⟩
    # = exp(-η^2 / 2) η^Δn √(γ(n₋ + 1) / γ(n₊ + 1)) L^Δn_n₋(η^2)
    # = exp(-η^2 / 2 + Δn log(η) + lγ(n₋ + 1) / 2 - lγ(n₊ + 1) / 2) L^Δn_n₋(η^2)
    η²::T = η * η
    @fastmath if n1 == n2
        lpre = η² * T(-0.5)
        lag = Utils.genlaguerre(n1, 0, η²)
    else
        n₋ = min(n1, n2)
        n₊ = max(n1, n2)
        Δn = abs(n1 - n2)
        lpre = (-η² + lgamma(T(n₋ + 1)) - lgamma(T(n₊ + 1))) / 2 + log(η) * Δn
        lag = Utils.genlaguerre(n₋, Δn, η²)
    end
    return @fastmath lag * exp(lpre)
end
sideband(n1::Integer, n2::Integer, η) = sideband(n1, n2, float(η))

@noinline function resize_caches_thread(tid, Ωcache, pcache)
    old_len = size(Ωcache, 1)
    resize!(Ωcache, tid)
    resize!(pcache, tid)
    for i in (old_len + 1):tid
        Ωcache[i] = Vector{Float64}(undef, 1024)
        pcache[i] = Vector{Float64}(undef, 1024)
    end
end

@generated function thermal_sideband(nbar::NTuple{N}, t, η::NTuple{N},
                                     Δn::NTuple{N}) where N
    factor = ceil(Int, log(N * 1e3))
    init = quote
        nmax = ($((:(ceil(Int, nbar[$i] * $factor)) for i in 1:N)...),)
        nbarp1 = ($((:(1 / (nbar[$i] + 1)) for i in 1:N)...),)
        α = ($((:(nbar[$i] * nbarp1[$i]) for i in 1:N)...),)
        s = 0.0
    end
    if N >= 2
        push!(init.args, :(tid = Threads.threadid()))
    end

    @gensym p
    @gensym Ω
    body = :(s += $p * sin($Ω)^2)

    for i in N:-1:2
        @gensym j
        @gensym p0
        @gensym Ωlocal
        @gensym plocal
        Ωcache = [Vector{Float64}(undef, 1024)]
        pcache = [Vector{Float64}(undef, 1024)]
        setup = quote
            if tid < size($Ωcache, 1)
                resize_caches_thread(tid, $Ωcache, $pcache)
            end
            @inbounds begin
                $Ωlocal = $Ωcache[tid]
                $plocal = $pcache[tid]
            end
            if size($Ωlocal, 1) < nmax[$i] + 1
                resize!($Ωlocal, nmax[$i] + 1)
                resize!($plocal, nmax[$i] + 1)
            end
            local $p0 = nbarp1[$i]
            @inbounds for $j in 0:nmax[$i]
                local $j
                $Ωlocal[$j + 1] = sideband($j, $j + Δn[$i], η[$i])
                $plocal[$j + 1] = $p0
                $p0 = $p0 * α[$i]
            end
        end
        push!(init.args, setup)

        prev_p = p
        prev_Ω = Ω
        @gensym p
        @gensym Ω
        body = quote
            for $j = 0:nmax[$i]
                local $j
                $prev_p = $p * $plocal[$j + 1]
                $prev_Ω = $Ω * $Ωlocal[$j + 1]
                $body
            end
        end
    end

    @gensym j
    @gensym p0
    body = quote
        local $p0 = nbarp1[1]
        @inbounds @fastmath for $j = 0:nmax[1]
            local $j
            $Ω = sideband($j, $j + Δn[1], η[1]) * t
            $p = $p0
            $p0 = $p0 * α[1]
            $body
        end
    end
    return quote
        $init
        $body
        return s
    end
end
thermal_sideband(nbar, t, η, Δn) = thermal_sideband((nbar,), t, (η,), (Δn,))

end
