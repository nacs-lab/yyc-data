#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: binomial_estimate
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170921_004940.csv")

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
                    if counts[idx, 1] != 0
                        loads[pg] += 1
                        if counts[idx, 2] != 0
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

idxs, ratios, uncs = group_time(params, counts, (1, 2, 3), 812 * 8)

const prefix = joinpath(@__DIR__, "imgs", "data_20170921_004940_cool_time")

figure()
errorbar(idxs, ratios[:, 1], uncs[:, 1], fmt=".-", label="X")
errorbar(idxs, ratios[:, 2], uncs[:, 2], fmt=".-", label="Y")
errorbar(idxs, ratios[:, 3], uncs[:, 3], fmt=".-", label="Z")
grid()
ylim([0, 0.14])
legend()
xlabel("Initial sequence index")
ylabel("Total cooling sideband")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
