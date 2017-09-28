#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: binomial_estimate
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170927_193235.csv")

data = readdlm(iname_a, ',', Float64, skipstart=1)
const params = data[:, 1]
const counts = Int.(@view data[:, 2:end])

function group_time(params, counts, param_grps, idx_grp_sz)
    np = length(params)
    ng = np ÷ idx_grp_sz
    npg = length(param_grps)
    ratios = Matrix{Float64}(ng, npg)
    uncs = Matrix{Float64}(ng, npg)
    loads = Vector{Int}(npg)
    survives = Vector{Int}(npg)
    for g in 1:ng
        loads .= 0
        survives .= 0
        for i in 1:idx_grp_sz
            idx = i + (g - 1) * idx_grp_sz
            param = params[idx]
            for pg in 1:npg
                if param in param_grps[pg]
                    if counts[idx, 2] != 0
                        loads[pg] += 1
                        if counts[idx, 3] != 0
                            survives[pg] += 1
                        end
                    end
                    break
                end
            end
        end
        for pg in 1:npg
            ratios[g, pg], uncs[g, pg] = binomial_estimate(survives[pg], loads[pg])
        end
    end
    [(0:(ng - 1)) .* idx_grp_sz .+ 1;], ratios, uncs
end

group1 = [1:5, 6:10, 11:15, 16:16]
group2 = [idx .+ 16 for idx in group1]
group3 = [idx .+ 16 for idx in group2]
group4 = [idx .+ 16 for idx in group3]

idxs, ratios, uncs = group_time(params, counts, [group1; group2; group3; group4], 800 * 16)

const prefix = joinpath(@__DIR__, "imgs", "data_20170927_193235_cool_time")

figure()
errorbar(idxs, ratios[:, 1], uncs[:, 1], fmt=".-", label="X")
errorbar(idxs, ratios[:, 2], uncs[:, 2], fmt=".-", label="Y")
errorbar(idxs, ratios[:, 3], uncs[:, 3], fmt=".-", label="Z")
grid()
ylim([0, 0.20])
legend()
xlabel("Initial sequence index")
ylabel("Total cooling sideband")
title("0.8")
NaCsPlot.maybe_save("$(prefix)_0.8")

figure()
errorbar(idxs, ratios[:, 5], uncs[:, 5], fmt=".-", label="X")
errorbar(idxs, ratios[:, 6], uncs[:, 6], fmt=".-", label="Y")
errorbar(idxs, ratios[:, 7], uncs[:, 7], fmt=".-", label="Z")
grid()
ylim([0, 0.20])
legend()
xlabel("Initial sequence index")
ylabel("Total cooling sideband")
title("1.0")
NaCsPlot.maybe_save("$(prefix)_1.0")

figure()
errorbar(idxs, ratios[:, 9], uncs[:, 9], fmt=".-", label="X")
errorbar(idxs, ratios[:, 10], uncs[:, 10], fmt=".-", label="Y")
errorbar(idxs, ratios[:, 11], uncs[:, 11], fmt=".-", label="Z")
grid()
ylim([0, 0.24])
legend()
xlabel("Initial sequence index")
ylabel("Total cooling sideband")
title("1.25")
NaCsPlot.maybe_save("$(prefix)_1.25")

figure()
errorbar(idxs, ratios[:, 13], uncs[:, 13], fmt=".-", label="X")
errorbar(idxs, ratios[:, 14], uncs[:, 14], fmt=".-", label="Y")
errorbar(idxs, ratios[:, 15], uncs[:, 15], fmt=".-", label="Z")
grid()
ylim([0, 0.28])
legend()
xlabel("Initial sequence index")
ylabel("Total cooling sideband")
title("1.5625")
NaCsPlot.maybe_save("$(prefix)_1.5625")

figure()
errorbar(idxs, (ratios[:, 1] .+ ratios[:, 2] .+ ratios[:, 3]) ./ 3,
         sqrt.(uncs[:, 1].^2 .+ uncs[:, 2].^2 .+ uncs[:, 3].^2) ./ 3, fmt=".-",
         label="0.8")
errorbar(idxs, (ratios[:, 5] .+ ratios[:, 6] .+ ratios[:, 7]) ./ 3,
         sqrt.(uncs[:, 5].^2 .+ uncs[:, 6].^2 .+ uncs[:, 7].^2) ./ 3, fmt=".-",
         label="1.0")
errorbar(idxs, (ratios[:, 9] .+ ratios[:, 10] .+ ratios[:, 11]) ./ 3,
         sqrt.(uncs[:, 9].^2 .+ uncs[:, 10].^2 .+ uncs[:, 11].^2) ./ 3, fmt=".-",
         label="1.25")
errorbar(idxs, (ratios[:, 13] .+ ratios[:, 14] .+ ratios[:, 15]) ./ 3,
         sqrt.(uncs[:, 13].^2 .+ uncs[:, 14].^2 .+ uncs[:, 15].^2) ./ 3, fmt=".-",
         label="1.5625")
grid()
ylim([0, 0.21])
xlabel("Initial sequence index")
ylabel("Total cooling sideband")
title("Average XYZ")
legend()
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(idxs, ratios[:, 4], uncs[:, 4], fmt=".-", label="0.8")
errorbar(idxs, ratios[:, 8], uncs[:, 8], fmt=".-", label="1.0")
errorbar(idxs, ratios[:, 12], uncs[:, 12], fmt=".-", label="1.25")
errorbar(idxs, ratios[:, 16], uncs[:, 16], fmt=".-", label="1.5625")
grid()
ylim([0.6, 1])
xlabel("Initial sequence index")
ylabel("Total survival")
legend()
NaCsPlot.maybe_save("$(prefix)_survival")

NaCsPlot.maybe_show()
