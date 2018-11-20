#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181118_214643.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const spec_a = [298.1292, 298.130275, 298.13098,
                298.1292, 298.130275, 298.13098,
                298.1292, 298.130275, 298.13098,
                298.1292, 298.130275, 298.13098, 0.0]
const iname_b = joinpath(@__DIR__, "data", "data_20181118_144941.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const spec_b = [298.129, 298.12985, 298.1307,
                298.129, 298.12985, 298.1307,
                298.129, 298.12985, 298.1307,
                298.129, 298.12985, 298.1307, 0.0]

function select_times_a(selector, binsz)
    nparam = length(params_a)
    return [NaCsData.split_data(NaCsData.select_count(params_a, logicals_a, selector,
                                                      i->i in i0:(i0 + binsz - 1)),
                                spec_a) for i0 in 1:binsz:(nparam - binsz ÷ 2)]
end
function select_times_b(selector, binsz)
    nparam = length(params_b)
    return [NaCsData.split_data(NaCsData.select_count(params_b, logicals_b, selector,
                                                      i->i in i0:(i0 + binsz - 1)),
                                spec_b) for i0 in 1:binsz:(nparam - binsz ÷ 2)]
end

const data_nacs_a = select_times_a(NaCsData.select_single((1, 2,), (3, 4,)), 1209 * 4)
const data_nacs_b = select_times_b(NaCsData.select_single((1, 2,), (3, 4,)), 1209 * 4)

const prefix = joinpath(@__DIR__, "imgs", "data_20181118_214643_raman_fluctuation")

function get_time_img(datas)
    nd = length(datas)
    local xs, img, img_unc
    for i in 1:nd
        data = [datas[i];]
        params, ratios, uncs = NaCsData.get_values(data)
        perm = sortperm(params)
        if !@isdefined xs
            xs = params[perm]
            img = zeros(nd, length(params))
            img_unc = zeros(nd, length(params))
        end
        img[i, :] = ratios[perm, 2]
        img_unc[i, :] = uncs[perm, 2]
    end
    return xs, img, img_unc
end

xs, img, unc = get_time_img(data_nacs_a)
figure()
for i in 2:length(xs)
    errorbar(1:size(img, 1), img[:, i], unc[:, i], label="$(i - 1)", fmt="o-")
end
ylim([0.3, 0.8])
legend(ncol=3)
grid()
xlabel("Time index")
ylabel("Suvival")
NaCsPlot.maybe_save("$(prefix)_time1")

figure()
errorbar(1:size(img, 1), img[:, 4] - img[:, 2], sqrt.(unc[:, 2].^2 .+ unc[:, 4].^2), fmt="o-")
grid()
xlabel("Time index")
ylabel("Suvival difference")
NaCsPlot.maybe_save("$(prefix)_diff1")

xs, img, unc = get_time_img(data_nacs_b)
figure()
for i in 2:length(xs)
    errorbar(1:size(img, 1), img[:, i], unc[:, i], label="$(i - 1)", fmt="o-")
end
ylim([0.3, 0.8])
legend(ncol=3)
grid()
xlabel("Time index")
ylabel("Suvival")
NaCsPlot.maybe_save("$(prefix)_time0")

figure()
errorbar(1:size(img, 1), img[:, 4] - img[:, 2], sqrt.(unc[:, 2].^2 .+ unc[:, 4].^2), fmt="o-")
grid()
xlabel("Time index")
ylabel("Suvival difference")
NaCsPlot.maybe_save("$(prefix)_diff0")

NaCsPlot.maybe_show()
