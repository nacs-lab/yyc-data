#!/usr/bin/julia

const δ0s = -50e3:500:50e3
const outdir = ARGS[1]
const prange = parse.(Int, split(ARGS[2], ','))
@assert all(0 .<= prange .< 8)

@everywhere include("clock_shift.jl")

@sync @distributed for p in prange
    BLAS.set_num_threads(2)
    run(p, outdir, δ0s, f1=f_cs, f2=f_na, z1=0.4, z2=1.0, maxf=900e3)
end
