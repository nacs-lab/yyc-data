#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190124_234335.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const binszs = [1200]
const specs_na = [1:2]
function select_times(params, logicals, selector, binsz, spec)
    nparam = length(params)
    return [NaCsData.split_data(NaCsData.select_count(params, logicals, selector,
                                                      i->i in i0:(i0 + binsz - 1)),
                                spec) for i0 in 1:binsz:(nparam - binsz ÷ 2)]
end
select_datas(datas, selector, binszs, specs) =
    [select_times(data..., selector, binsz, spec)
     for (data, binsz, spec) in zip(datas, binszs, specs)]
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), binszs, specs_na)

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

const prefix = joinpath(@__DIR__, "imgs", "data_20190124_234335_depump_rt")

xs, img, unc = get_time_img(datas_na[1])

figure()
errorbar(1:size(img, 1), img[:, 1], unc[:, 1], label="2, 2", fmt="o-")
errorbar(1:size(img, 1), img[:, 2], unc[:, 2], label="2, 1", fmt="o-")
legend()
grid()
xlabel("Time index")
ylabel("Suvival")
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(1:size(img, 1), img[:, 1] ./ img[:, 2],
         sqrt.(unc[:, 1].^2 ./ img[:, 2].^2 .+ unc[:, 2].^2 .* img[:, 1].^2 ./ img[:, 2].^4),
         fmt="o-")
grid()
xlabel("Time index")
ylabel("Suvival")
NaCsPlot.maybe_save("$(prefix)_ratio")

NaCsPlot.maybe_show()
