#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181120_010424.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const spec_a = (0.0:0.0, [0.0, 1.0], 0.2:0.2:2.0)

function select_times(selector, binsz)
    nparam = length(params_a)
    datas = [NaCsData.split_data(NaCsData.select_count(params_a, logicals_a, selector,
                                                      i->i in i0:(i0 + binsz - 1)),
                                spec_a) for i0 in 1:binsz:(nparam - binsz ÷ 2)]
    return [[data[1]; data[3]] for data in datas], [data[2] for data in datas]
end

const data_nacs_a = select_times(NaCsData.select_single((1, 2,), (3, 4,)), 1230 * 4)

const prefix = joinpath(@__DIR__, "imgs", "data_20181120_010424_raman_time_rt")

function get_time_img(datas)
    nd = length(datas)
    local xs, img
    for i in 1:nd
        data = datas[i]
        params, ratios, uncs = NaCsData.get_values(data)
        perm = sortperm(params)
        if !@isdefined xs
            xs = params[perm]
            img = zeros(nd, length(params))
        end
        img[i, :] = ratios[perm, 2]
    end
    return xs, img
end

xs, img = get_time_img(data_nacs_a[1])
figure()
imshow(img, interpolation="gaussian", extent=[xs[1], xs[end], size(img, 1), 1];
       cmap="viridis", aspect="auto")
title("Raman Time")
xlabel("Time (ms)")
ylabel("Time index")
colorbar()
NaCsPlot.maybe_save("$(prefix)_time")

xs, img = get_time_img(data_nacs_a[2])
figure()
imshow(img, interpolation="gaussian", extent=[xs[1], xs[end], size(img, 1), 1];
       cmap="viridis", aspect="auto")
title("Reference points")
ylabel("Time index")
colorbar()
NaCsPlot.maybe_save("$(prefix)_ref")

NaCsPlot.maybe_show()
