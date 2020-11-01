#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using LinearAlgebra

const freqs10 = [4.0 1.0 103.630 0.05
                 4.0 1.0 105.700 0.05
                 4.0 1.0 107.700 0.1
                 6.0 1.0 105.320 0.06
                 6.0 1.0 107.180 0.06
                 6.0 1.0 108.164 0.06
                 7.3 1.0 106.200 0.06
                 7.3 1.0 108.000 0.06
                 7.3 1.0 109.000 0.06
                 8.6 1.0 106.740 0.06
                 8.6 1.0 108.400 0.06
                 8.6 1.0 110.260 0.10
                 10  1.0 107.340 0.06
                 10  1.0 109.020 0.06
                 10  1.0 111.720 0.10
                 12  1.0 107.800 0.10
                 12  1.0 109.864 0.06
                 12  1.0 113.920 0.10]
const freqs06 = [2 0.6 101.480 0.05
                 2 0.6 102.816 0.05
                 2 0.6 104.500 0.12
                 3 0.6 102.360 0.05
                 3 0.6 103.600 0.05
                 3 0.6 104.400 0.05
                 4 0.6 103.050 0.05
                 4 0.6 104.200 0.05
                 4 0.6 104.860 0.05
                 5 0.6 103.560 0.06
                 5 0.6 104.688 0.06
                 5 0.6 105.840 0.06
                 6 0.6 103.924 0.06
                 6 0.6 105.088 0.06
                 6 0.6 106.960 0.06
                 7 0.6 104.160 0.12
                 7 0.6 105.600 0.06
                 7 0.6 107.800 0.12
                 8 0.6 104.320 0.12
                 8 0.6 106.060 0.12
                 8 0.6 108.800 0.12]

const freqs = [freqs10; freqs06]

# Clebsh-Gordan coeeficients
#           ⟨NF=8, m_NF=6| ⟨NF=7, m_NF=6| ⟨NF=6, m_NF=6|
const ce = [sqrt(1 / 20)  -sqrt(9 / 28)   sqrt(22 / 35) # ⟨N=2 m_N=0; F=6 m_F=6|
            sqrt(2 / 5)   -sqrt(2 / 7)   -sqrt(11 / 35) # ⟨N=2 m_N=1; F=6 m_F=5|
            sqrt(11 / 20)  sqrt(11 / 28)  sqrt(2 / 35)] # ⟨N=2 m_N=2; F=6 m_F=4|
const M1 = ce[:, 1] * ce[:, 1]'
const M2 = ce[:, 2] * ce[:, 2]'
const M3 = ce[:, 3] * ce[:, 3]'

# function elevels(power, B, ps)
#     M = M1 .* -(ps[2] + ps[3]) .+ M2 .* ps[2] .+ M3 .* ps[3]
#     μB = B * ps[7]
#     M[1, 1] += ps[1] + power * ps[4]
#     M[2, 2] += ps[1] + power * ps[5] + μB
#     M[3, 3] += ps[1] + power * ps[6] + μB * 2
#     return eigvals(M)
# end
function elevels(power, B, ps)
    M = M1 .* -(ps[2] + ps[3]) .+ M2 .* ps[2] .+ M3 .* ps[3]
    # μB = B * ps[7]
    M[1, 1] += ps[1] + power * ps[4] + B * ps[7]
    M[2, 2] += ps[1] + power * ps[5] + B * ps[8]
    M[3, 3] += ps[1] + power * ps[6] + B * ps[9]
    return eigvals(M)
end

function model(x, ps)
    cache = Dict{Tuple{Float64,Float64},Vector{Float64}}()
    Es = Float64[]
    for idx in x
        power = freqs[idx, 1]
        B = freqs[idx, 2]
        ls = get!(cache, (power, B)) do
            elevels(power, B, ps)
        end
        push!(Es, ls[(idx - 1) % 3 + 1])
    end
    return Es
end
fit = fit_data(model, 1:size(freqs, 1), freqs[:, 3], [99.2, 0.0, 1.26,
                                                      1.15, 0.52, 0.14,
                                                      0, 4.16, 8.32], plotx=false)
@show fit.uncs
@show elevels(0, 0, fit.param)

function eval_model_power(powers, B, ps)
    E1s = Float64[]
    E2s = Float64[]
    E3s = Float64[]
    for power in powers
        E = elevels(power, B, ps)
        append!(E1s, E[1])
        append!(E2s, E[2])
        append!(E3s, E[3])
    end
    return powers, E1s, E2s, E3s
end

function eval_model_B(power, Bs, ps)
    E1s = Float64[]
    E2s = Float64[]
    E3s = Float64[]
    for B in Bs
        E = elevels(power, B, ps)
        append!(E1s, E[1])
        append!(E2s, E[2])
        append!(E3s, E[3])
    end
    return Bs, E1s, E2s, E3s
end

const prefix = joinpath(@__DIR__, "imgs", "fit_20181103_raman_n2")

plotx, plote1, plote2, plote3 = eval_model_power(linspace(3.5, 13, 1000), 1.0, fit.param)

figure()
plot(plotx, plote1, "C1")
plot(plotx, plote2, "C1")
plot(plotx, plote3, "C1")
errorbar(freqs10[:, 1], freqs10[:, 3], freqs10[:, 4], fmt="C0o")
grid()
title("Full B field")
xlabel("Power (mW)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_b10")

plotx, plote1, plote2, plote3 = eval_model_power(linspace(1.5, 10, 1000), 0.6, fit.param)

figure()
plot(plotx, plote1, "C1")
plot(plotx, plote2, "C1")
plot(plotx, plote3, "C1")
errorbar(freqs06[:, 1], freqs06[:, 3], freqs06[:, 4], fmt="C0o")
grid()
title("60% B field")
xlabel("Power (mW)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_b06")

plotx, plote1, plote2, plote3 = eval_model_B(0, linspace(0, 1, 1001), fit.param)

figure()
plot(plotx .* 8.85, plote1, "C1")
plot(plotx .* 8.85, plote2, "C1")
plot(plotx .* 8.85, plote3, "C1")
grid()
title("0 Power")
xlabel("B field (G)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_p0")

NaCsPlot.maybe_show()
