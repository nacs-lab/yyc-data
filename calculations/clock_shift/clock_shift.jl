#

include("utils.jl")

if VERSION >= v"0.7.0"
    @eval using LinearAlgebra
end

# Cs traping frequencies 20, 130, 140
# Na traping frequencies 84% that of Cs

const f_cs = (20e3, 130e3, 140e3)
const f_na = f_cs .* 0.84

function iterate_ax(cb, f, maxf)
    for i in 0:typemax(Int)
        δ = maxf - f * i
        if δ < 0
            return
        end
        cb(i, δ)
    end
end

function iterate_3ax(cb, fs, maxf)
    iterate_ax(fs[1], maxf) do i1, maxf1
        iterate_ax(fs[2], maxf1) do i2, maxf2
            iterate_ax(fs[3], maxf2) do i3, maxf3
                cb((i1, i2, i3), maxf3)
            end
        end
    end
end

function iterate_2atoms(cb, fs1, fs2, maxf)
    iterate_3ax(fs1, maxf) do i1, maxf1
        iterate_3ax(fs2, maxf1) do i2, δ
            cb((i1..., i2...), δ)
        end
    end
end

function count_states(iterf, args...)
    c = 0
    iterf(args...) do _1, _2
        c += 1
    end
    return c
end

function collect_states(iterf, args...)
    local res
    iterf(args...) do i, _2
        if @isdefined(res)
            push!(res, i)
        else
            res = [i]
        end
    end
    return res
end

struct WFOverlapCache
    z1::Float64
    z2::Float64
    cache::Dict{NTuple{4,Int},Float64}
    WFOverlapCache(z1, z2) = new(z1, z2, Dict{NTuple{4,Int},Float64}())
end

function (cache::WFOverlapCache)(n1::Integer, n2::Integer, m1::Integer, m2::Integer)
    if (n1 + n2 + m1 + m2) % 2 != 0
        return 0.0
    end
    ks = keys(cache.cache)
    if (n1, n2, m1, m2) in ks
        return cache.cache[(n1, n2, m1, m2)]
    elseif (n2, n1, m1, m2) in ks
        return cache.cache[(n2, n1, m1, m2)]
    elseif (n1, n2, m2, m1) in ks
        return cache.cache[(n1, n2, m2, m1)]
    elseif (n2, n1, m2, m1) in ks
        return cache.cache[(n2, n1, m2, m1)]
    end
    v = wavefunction_overlap(n1, n2, m1, m2, cache.z1, cache.z2)
    cache.cache[(n1, n2, m1, m2)] = v
    return v
end

function populate_matrix(fs1, fs2, maxf, z1, z2, δ0)
    states = collect_states(iterate_2atoms, fs1, fs2, maxf)
    n = length(states)
    fs = (fs1..., fs2...)
    res = @static VERSION >= v"0.7.0" ? Matrix{Float64}(undef, n, n) : Matrix{Float64}(n, n)
    cache = WFOverlapCache(z1, z2)
    δ0 = δ0 / cache(0, 0, 0, 0)^3
    for i in 1:n
        state1 = states[i]
        for j in i:n
            if (i == j)
                v = sum(fs .* state1)
            else
                v = 0.0
            end
            state2 = states[j]
            v += (cache(state1[1], state2[1], state1[4], state2[4]) *
                  cache(state1[2], state2[2], state1[5], state2[5]) *
                  cache(state1[3], state2[3], state1[6], state2[6])) * δ0
            res[i, j] = v
        end
    end
    return Symmetric(res)
end
