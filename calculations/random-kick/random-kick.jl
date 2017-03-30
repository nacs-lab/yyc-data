#!/usr/bin/julia

@inline promote_float(args...) = float.(promote(args...))

function propagate(tlen::T, F::T, δ::T, decay::T) where {T<:AbstractFloat}
    A = zero(Complex{T})
    t = float(zero(tlen))
    if 0 < decay < 1
        nmax = ceil(Int, log(rand()) / log(1 - decay))
    else
        nmax = typemax(Int)
    end
    nkick = 0
    while true
        r = T(rand() * 2 - 1)
        r == 0 && break
        ra = abs(r)
        δt = @fastmath -log(ra) / F
        if t + δt < tlen
            A *= @fastmath cis(δt)
            if r > 0
                A += 2δ
            end
            t += δt
            nkick += 1
            if nkick >= nmax
                break
            end
        else
            break
        end
    end
    return abs2(A)
end

function sample(N::Integer, tlen::T, F::T, δ::T, decay::T) where {T<:AbstractFloat}
    ΣA = zero(T)
    ΣA2 = zero(T)
    for i in 1:N
        A = propagate(tlen, F, δ, decay)
        ΣA += A
        ΣA2 += A^2
    end
    avgA = ΣA / N
    avgA2 = ΣA2 / N
    return avgA, sqrt((avgA2 - avgA^2) / (N - 1))
end

function samples(N, ts::AbstractArray{T}, F::T, δ::T, decay::T) where {T<:AbstractFloat}
    nt = length(ts)
    rs = Vector{T}(nt)
    ss = Vector{T}(nt)
    @inbounds for i in 1:nt
        res = sample(N, ts[i], F, δ, decay)
        rs[i] = res[1]
        ss[i] = res[2]
    end
    return rs, ss
end

samples(N, ts::AbstractArray{T}, F, δ, decay=0) where {T<:AbstractFloat} =
    samples(N, ts, T(F), T(δ), T(decay))
function samples(N, ts::AbstractArray, _F, _δ, _decay=0)
    F, δ, decay = promote_float(_F, _δ, _decay)
    T = typeof(F)
    return samples(N, T.(ts), F, δ, decay)
end

using PyPlot
ts = 0:0.1f0:10f0
decay = 0.2f0
errorbar(ts, samples(1000000, ts / 0.01f0, 0.01f0, 1f0, decay)..., label="0.01")
errorbar(ts, samples(1000000, ts / 0.423f0, 0.423f0, 1f0, decay)..., label="0.423")
errorbar(ts, samples(1000000, ts / 0.146f0, 0.146f0, 1f0, decay)..., label="0.146")
errorbar(ts, samples(1000000, ts / 0.23f0, 0.23f0, 1f0, decay)..., label="0.23")
grid()
legend()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
show()
