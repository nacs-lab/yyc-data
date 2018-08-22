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
    cutoff::Float64
    cache::Dict{NTuple{4,Int},Float64}
    WFOverlapCache(z1, z2, cutoff) = new(z1, z2, abs(cutoff), Dict{NTuple{4,Int},Float64}())
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
    if v > cache.cutoff
        v = cache.cutoff
    elseif v < -cache.cutoff
        v = -cache.cutoff
    end
    cache.cache[(n1, n2, m1, m2)] = v
    return v
end

# p0 is the parity of each axis
# The interaction does not mix parity on each axis so we can compute those separately
# This reduce the number of states by 8. Since the memory scales with n^2 and time with n^3
# this is a big win.
function coupled_2atoms(fs1, fs2, maxf, p0; maxtotaln=-1, maxns=())
    states = NTuple{6,Int}[]
    if maxtotaln < 0
        maxtotaln = typemax(Int)
    else
        maxtotaln = Int(maxtotaln)
    end
    _maxns = ntuple(i->i > length(maxns) ? typemax(Int) : Int(maxns[i]), 6)
    iterate_2atoms(fs1, fs2, maxf) do i, _2
        if ((i[1] + i[4] - p0[1]) % 2 != 0 || (i[2] + i[5] - p0[2]) % 2 != 0 ||
            (i[3] + i[6] - p0[3]) % 2 != 0)
            return
        end
        if sum(i) > maxtotaln || any(i .> _maxns)
            return
        end
        push!(states, i)
    end
    return states
end

struct H2Atoms
    states::Vector{NTuple{6,Int}}
    es::Vector{Float64}
    inter::Matrix{Float64}
end

function H2Atoms(fs1, fs2, maxf, z1, z2, p0; cutoff=Inf, maxtotaln=-1, maxns=())
    states = coupled_2atoms(fs1, fs2, maxf, p0, maxtotaln=maxtotaln, maxns=maxns)
    n = length(states)
    fs = (fs1..., fs2...)
    inter = @static VERSION >= v"0.7.0" ? Matrix{Float64}(undef, n, n) : Matrix{Float64}(n, n)
    es = @static VERSION >= v"0.7.0" ? Vector{Float64}(undef, n) : Vector{Float64}(n)
    cache = WFOverlapCache(z1, z2, cutoff)
    for i in 1:n
        state1 = states[i]
        es[i] = sum(fs .* state1)
        for j in i:n
            state2 = states[j]
            v = (cache(state1[1], state2[1], state1[4], state2[4]) *
                 cache(state1[2], state2[2], state1[5], state2[5]) *
                 cache(state1[3], state2[3], state1[6], state2[6]))
            inter[i, j] = v
        end
    end
    return H2Atoms(states, es, inter)
end

function getH(h0::H2Atoms, δ0)
    δ0 = δ0 / h0.inter[1, 1]
    n = length(h0.es)
    H = @static VERSION >= v"0.7.0" ? Matrix{Float64}(undef, n, n) : Matrix{Float64}(n, n)
    @inbounds for i in 1:n
        @simd for j in i:n
            H[i, j] = h0.inter[i, j] * δ0
        end
        H[i, i] += h0.es[i]
    end
    return Symmetric(H)
end
