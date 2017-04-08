#!/usr/bin/julia

import GSL

function coupling_LS(L1x2::Integer, L2x2::Integer, J1x2::Integer, J2x2::Integer, Ix2::Integer,
                     F1x2::Integer, F2x2::Integer, mF1x2::Integer, mF2x2::Integer)::Float64
    if abs(L1x2 - L2x2) != 2
        return 0
    elseif abs(J1x2 - L1x2) != 1 || abs(J2x2 - L2x2) != 1
        throw(ArgumentError("Illegal J"))
    elseif !(abs(J1x2 - Ix2) <= F1x2 <= J1x2 + Ix2 && abs(J2x2 - Ix2) <= F2x2 <= J2x2 + Ix2)
        throw(ArgumentError("Illegal F"))
    elseif (F1x2 + J1x2 + Ix2) % 2 != 0 || (F2x2 + J2x2 + Ix2) % 2 != 0
        throw(ArgumentError("Illegal F"))
    elseif !((F1x2 - F2x2) in (-2, 0, 2) && (mF1x2 - mF2x2) in (-2, 0, 2))
        return 0
    end
    qx2 = mF1x2 - mF2x2
    m1powx2 = F2x2 + mF1x2 - 2
    factor = sqrt(F1x2 + 1) * GSL.sf_coupling_3j(F2x2, 2, F1x2, mF2x2, qx2, -mF1x2)
    m1powx2 += F2x2 + J1x2 + 2 + Ix2
    factor *= sqrt((F2x2 + 1) * (J1x2 + 1)) * GSL.sf_coupling_6j(J1x2, J2x2, 2, F2x2, F1x2, Ix2)
    m1powx2 += J2x2 + L1x2 + 2 + 1
    factor *= sqrt((J2x2 + 1) * (L1x2 + 1)) * GSL.sf_coupling_6j(L1x2, L2x2, 2, J2x2, J1x2, 1)
    if m1powx2 % 4 == 2
        factor = -factor
    end
    return factor
end

coupling_D(D2::Bool, Ix2, F1x2, F2x2, mF1x2, mF2x2) =
    coupling_LS(0, 2, 1, D2 ? 3 : 1, Ix2, F1x2, F2x2, mF1x2, mF2x2)

function na_scattering_d2(F1::Integer, F2::Integer, mF1::Integer, mF2::Integer, pol::Integer)
    if !(F1 in (1, 2)) || !(F2 in (1, 2))
        throw(ArgumentError("Invalid F for Na"))
    elseif !(pol in (-1, 0, 1))
        throw(ArgumentError("Invalid polarization"))
    end
    mF′ = mF1 + pol

    F1x2 = F1 * 2
    mF1x2 = mF1 * 2
    F2x2 = F2 * 2
    mF2x2 = mF2 * 2
    mF′x2 = mF′ * 2

    cpl = 0.0
    for F′ in (F1 - 1):(F1 + 1)
        if abs(mF′) > F′
            continue
        end
        if abs(F′ - F2) > 1
            continue
        end
        F′x2 = F′ * 2
        cpl += (coupling_D(true, 3, F1x2, F′x2, mF1x2, mF′x2) *
                coupling_D(true, 3, F2x2, F′x2, mF2x2, mF′x2))
    end
    return cpl
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
const rates_f2_counterop = gen_rates((0, 0, 1), -25.0e9 - 1.77e9)

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
