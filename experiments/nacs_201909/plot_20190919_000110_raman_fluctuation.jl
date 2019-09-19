#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190919_000110.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
# const maxcnts = [typemax(Int)]
const specs = [731.7 .+ [-0.2, 0.2]]
select_data(data, selector, maxcnt, spec) =
    NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
function select_time(data, selector, binsz, spec)
    nparam = length(data[1])
    return [NaCsData.split_data(NaCsData.select_count(data..., selector,
                                                      i->i in i0:(i0 + binsz - 1)),
                                spec) for i0 in 1:binsz:(nparam - binsz ÷ 2)]
end

select_times(datas, selector, binsz, specs) =
    [select_time(data, selector, binsz, spec) for (data, spec) in zip(datas, specs)]

const datas_nacs = select_times(datas, NaCsData.select_single((1, 2,), (3, 4,)),
                                1209 * 6, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20190919_000110_raman_fluctuation")

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

xs, img, unc = get_time_img(datas_nacs[1])
figure()
errorbar(1:size(img, 1), img[:, 1], unc[:, 1], label="Left", fmt="C0o-")
errorbar(1:size(img, 1), img[:, 2], unc[:, 2], label="Right", fmt="C1o-")
legend(ncol=2, fontsize="small")
grid()
xlabel("Time index")
ylabel("Suvival")
title("Left vs Right")
NaCsPlot.maybe_save("$(prefix)_lr")

figure()
errorbar(1:size(img, 1), img[:, 2] - img[:, 1], sqrt.(unc[:, 2].^2 .+ unc[:, 1].^2), fmt="o-")
grid()
xlabel("Time index")
ylabel("Suvival difference")
title("Difference")
NaCsPlot.maybe_save("$(prefix)_diff")

NaCsPlot.maybe_show()
