#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci
using MAT

const iname_a = joinpath(@__DIR__, "data", "data_20180108_181906.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const counts_p, counts_v, counts_u = matopen(iname_a) do f
    cs = read(f, "Counts")::Array{Float64,3}
    nimgs = size(cs, 1)
    nsites = size(cs, 2)
    data_dict = Dict{Float64,Tuple{Matrix{Float64},Matrix{Float64},typeof(Ref(0))}}()
    for seq in 1:length(params_a)
        param = params_a[seq]
        count = @view cs[:, :, seq]
        if haskey(data_dict, param)
            item = data_dict[param]
        else
            item = data_dict[param] = (zeros(nimgs, nsites), zeros(nimgs, nsites), Ref(0))
        end
        item[1] .+= count
        item[2] .+= count.^2
        item[3][] += 1
    end
    params = sort(collect(keys(data_dict)))
    len = length(params)
    avgs = zeros(nimgs, nsites, len)
    uncs = zeros(nimgs, nsites, len)
    for i in 1:len
        p = data_dict[params[i]]
        nrep = p[3][]
        for j in 1:nsites
            for k in 1:nimgs
                avg = p[1][k, j] / nrep
                avg2 = p[2][k, j] / nrep
                avgs[k, j, i] = avg
                uncs[k, j, i] = sqrt(avg2 - avg^2) / sqrt(nrep - 1)
            end
        end
    end
    params, avgs, uncs
end

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

const data = NaCsData.select_count(params_a, logicals_a, selector)

const prefix = joinpath(@__DIR__, "imgs", "data_20180108")

figure()
NaCsPlot.plot_loading_data(data, label="Loading")
NaCsPlot.plot_survival_data(data, label="Survival")
grid()
ylim([0, 1])
legend()
title("Probability")
xlabel("Raising edge (deg)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_prob")


figure()
errorbar(counts_p, counts_v[1, 1, :], counts_u[1, 1, :], label="First")
errorbar(counts_p, counts_v[3, 1, :], counts_u[3, 1, :], label="Second")
grid()
legend()
title("Counts")
xlabel("Raising edge (deg)")
ylabel("Signal")
NaCsPlot.maybe_save("$(prefix)_count")

NaCsPlot.maybe_show()
