#!/usr/bin/julia

include("clock_shift.jl")

const δ0s = -50e3:100:50e3

const outdir = ARGS[1]
const prange = parse.(Int, split(ARGS[2], ','))
@assert all(0 .<= prange .< 8)

function dump_sys(io, h::H2Atoms, p::NTuple{3,Int})
    @assert length(h.states) == length(h.es) == size(h.inter, 1) == size(h.inter, 2)
    @assert Int === Int64
    write(io, length(h.states))
    write(io, Ref(p))
    write(io, h.states)
    write(io, h.es)
    write(io, h.inter)
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

function run(p)
    p0 = ((p & 4) >> 2, (p & 2) >> 1, (p & 1) >> 0)
    @show p0
    h = H2Atoms(f_na, f_cs, 1000e3, 1, 0.4, p0,
                cutoff=1, maxtotaln=40, maxns=(29, 14, 14, 29, 14, 14))
    dir = joinpath(outdir, "$p")
    mkpath(dir)
    open(joinpath(dir, "sys.bin"), "w") do io
        dump_sys(io, h, p0)
    end
    n = length(h.states)
    H = Symmetric(@static VERSION >= v"0.7.0" ?
                  Matrix{Float64}(undef, n, n) :
                  Matrix{Float64}(n, n))
    for i in 1:length(δ0s)
        δ0 = δ0s[i]
        println("$i/$(length(δ0s)): δ0 = $δ0")
        @time res = eigen(getH(h, δ0, H), 1:50)
        open(joinpath(dir, "res_$δ0.bin"), "w") do io
            dump_res(io, δ0, res.values, res.vectors)
        end
    end
end

for p in prange
    run(p)
end
