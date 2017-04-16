#!/usr/bin/julia

using NaCsCalc.Atomic: coupling_D_raman

function na_scattering_d2(F1::Integer, F2::Integer, mF1::Integer, mF2::Integer, pol::Integer)
    return coupling_D_raman(true, 3, F1 * 2, F2 * 2, mF1 * 2, mF2 * 2, pol)
end

function na_scattering_d2_rate(F1::Integer, F2::Integer, mF1::Integer, mF2::Integer,
                               pol::NTuple{3,Real}, δc::Real)
    if F1 == 1
        δ = δc + 664.360e6
    else
        δ = δc - 1.107266e9
    end
    rate = 0.0
    for i in 1:3
        p = pol[i]
        if p == 0
            continue
        end
        rate += na_scattering_d2(F1, F2, mF1, mF2, i - 2)^2 * p
    end
    return rate * (61.542e6 / δ)^2
end

const states = Tuple{Int,Int}[]
for i in -2:2
    push!(states, (2, i))
end
for i in -1:1
    push!(states, (1, i))
end
const nstates = length(states)

function gen_rates(pol, δ)
    rates = Matrix{Float64}(nstates, nstates)
    @inbounds for i in 1:nstates
        state1 = states[i]
        for j in 1:nstates
            state2 = states[j]
            rates[j, i] = na_scattering_d2_rate(state1[1], state2[1],
                                                state1[2], state2[2], pol, δ)
        end
    end
    return rates
end

const rates_f1_up = gen_rates((0.25, 0.5, 0.25), -25.0e9)
const rates_f1_down = rates_f1_up
const rates_f1_coprop = gen_rates((0.5, 0.0, 0.5), -25.0e9)
const rates_f2_coprop = gen_rates((0.25, 0.5, 0.25), -25.0e9 - 1.77e9)
const rates_f2_counterop = gen_rates((0.1, 0, 0.9), -25.0e9 - 1.77e9)

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
    res = expm(A * t) * init
    return res[6] + res[7] + res[8]
end

function gen_model(rates, init)
    A = rates_to_A(rates)
    t->propagate_f1(A, init, t)
end
