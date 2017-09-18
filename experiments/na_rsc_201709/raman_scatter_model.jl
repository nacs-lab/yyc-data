#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
using NaCsCalc.Atomic: all_scatter_D

const δf1 = -76.82e9
const δf2 = -76.82e9 - 1.77e9
const rlof_f1 = (61.542e6 / (δf1 - 1.107266e9))^2
const rlof_f2 = (61.542e6 / (δf2 - 1.107266e9))^2
const rhif_f1 = (61.542e6 / (δf1 + 664.360e6))^2
const rhif_f2 = (61.542e6 / (δf2 + 664.360e6))^2

const rates_f1_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f1_down = rates_f1_up
const rates_f1_diagonal = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f1, rlof_f1)
const rates_f2_diagonal = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f2, rlof_f2)
const rates_f2_counterop = all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_f2, rlof_f2)

function rates_to_A(rates)
    nx, ny = size(rates)
    A = Matrix{Float64}(nx, ny)
    @inbounds for i in 1:nx
        s = 0.0
        for j in 1:ny
            r = rates[j, i]
            A[j, i] = r
            s += r
        end
        A[i, i] -= s
    end
    return A
end

function propagate_f1(A, init, t)
    res = exp(A .* t) * init # expm!!!
    return res[6] + res[7] + res[8]
end

function gen_model(rates, init)
    A = rates_to_A(rates)
    t->propagate_f1(A, init, t)
end
