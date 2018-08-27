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
    δscale::Float64
end

function H2Atoms(fs1, fs2, maxf, z1, z2, p0; cutoff=Inf, maxtotaln=-1, maxns=())
    states = coupled_2atoms(fs1, fs2, maxf, p0, maxtotaln=maxtotaln, maxns=maxns)
    n = length(states)
    fs = (fs1..., fs2...)
    inter = uninit_ary(Matrix{Float64}, n, n)
    es = uninit_ary(Vector{Float64}, n)
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
    return H2Atoms(states, es, inter, wavefunction_overlap(0, 0, 0, 0, z1, z2)^3)
end

function getH(h0::H2Atoms, δ0, out)
    δ0 = δ0 / h0.δscale
    n = length(h0.es)
    H = parent(out)
    @inbounds for i in 1:n
        @simd for j in i:n
            H[i, j] = h0.inter[i, j] * δ0
        end
        H[i, i] += h0.es[i]
    end
    return out
end

function getH(h0::H2Atoms, δ0)
    n = length(h0.es)
    out = Symmetric(uninit_ary(Matrix{Float64}, n, n))
    return getH(h0, δ0, out)
end

function dump_sys(io, h::H2Atoms, p::NTuple{3,Int})
    @assert length(h.states) == length(h.es) == size(h.inter, 1) == size(h.inter, 2)
    @assert Int === Int64
    write(io, length(h.states))
    write(io, Ref(p))
    write(io, h.states)
    write(io, h.es)
    write(io, h.inter)
    write(io, h.δscale)
    return
end

function dump_res(io, δ0::Float64, vals::Vector{Float64}, vecs::Matrix{Float64})
    @assert length(vals) == size(vecs, 2)
    write(io, size(vecs, 1))
    write(io, size(vecs, 2))
    write(io, δ0)
    write(io, vals)
    write(io, vecs)
    return
end

function run(p, outdir, δ0s; f1=f_na, f2=f_cs, z1=1.0, z2=0.4, maxf=1000e3)
    p0 = ((p & 4) >> 2, (p & 2) >> 1, (p & 1) >> 0)
    @show p0
    h = H2Atoms(f1, f2, maxf, z1, z2, p0,
                cutoff=1, maxtotaln=40, maxns=(29, 14, 14, 29, 14, 14))
    dir = joinpath(outdir, "$p")
    mkpath(dir)
    open(joinpath(dir, "sys.bin"), "w") do io
        dump_sys(io, h, p0)
    end
    n = length(h.states)
    H = Symmetric(uninit_ary(Matrix{Float64}, n, n))
    for i in 1:length(δ0s)
        δ0 = δ0s[i]
        println("$i/$(length(δ0s)): δ0 = $δ0")
        @time res = eigen(getH(h, δ0, H), 1:50)
        open(joinpath(dir, "res_$δ0.bin"), "w") do io
            dump_res(io, δ0, res.values, res.vectors)
        end
    end
end

function load_sys(io)
    @assert Int === Int64
    n = read!(io, Ref{Int}())[]
    p = read!(io, Ref{NTuple{3,Int}}())[]
    states = uninit_ary(Vector{NTuple{6,Int}}, n)
    es = uninit_ary(Vector{Float64}, n)
    inter = uninit_ary(Matrix{Float64}, n, n)
    read!(io, states)
    read!(io, es)
    read!(io, inter)
    return p, H2Atoms(states, es, inter, 1) # Forgot to save δscale...
end

load_sys(fname::AbstractString) = open(load_sys, fname)

struct Res
    δ0::Float64
    vals::Vector{Float64}
    vecs::Matrix{Float64}
end

function load_res(io)
    ndim = read!(io, Ref{Int}())[]
    nres = read!(io, Ref{Int}())[]
    δ0 = read!(io, Ref{Float64}())[]
    vals = uninit_ary(Vector{Float64}, nres)
    vecs = uninit_ary(Matrix{Float64}, ndim, nres)
    read!(io, vals)
    read!(io, vecs)
    return Res(δ0, vals, vecs)
end
load_res(fname::AbstractString) = open(load_res, fname)

struct ResGroup
    p::NTuple{3,Int}
    h::H2Atoms
    res::Vector{Res}
    r0::Res
end

function load_dir(datadir)
    p, h = load_sys(joinpath(datadir, "sys.bin"))
    res = Res[]
    local r0
    for d in readdir(datadir)
        d == "sys.bin" && continue
        r = load_res(joinpath(datadir, d))
        push!(res, r)
        if r.δ0 == 0
            r0 = r
        end
    end
    sort!(res, by=r->r.δ0)
    @assert @isdefined(r0)
    return ResGroup(p, h, res, r0)
end

get_δ_it(rg::ResGroup) = (r.δ0 for r in rg.res)
get_energy_it(rg::ResGroup, i; base=0) = (r.vals[i] - base for r in rg.res)
get_overlap_it(rg::ResGroup, i, vec) = (abs(dot(r.vecs[:, i], vec)) for r in rg.res)

get_δ(rg::ResGroup) = collect(get_δ_it(rg))
get_energy(rg::ResGroup, i; base=0) = collect(get_energy_it(rg, i; base=base))
get_overlap(rg::ResGroup, i, vec) = collect(get_overlap_it(rg, i, vec))

function get_state(rg::ResGroup, ns::NTuple{6,Int})
    states = rg.h.states
    n = length(states)
    v = zeros(Float64, n)
    for i in 1:n
        if states[i] == ns
            v[i] = 1
            return v
        end
    end
    throw(ArgumentError("State $ns not found."))
end

# Max overlap between
max_overlap(rg::ResGroup, i, vec) = maximum(get_overlap_it(rg, i, vec))

filter_overlap(rg::ResGroup, vec, minovrlap) =
    (i for i in 1:length(rg.r0.vals) if max_overlap(rg, i, vec) >= minovrlap)
