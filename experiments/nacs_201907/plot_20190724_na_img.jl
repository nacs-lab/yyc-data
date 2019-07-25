#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const det = linspace(-40, -5, 8)
const amp = linspace(0.1, 0.8, 15)

const inames = ["data_20190724_182704.mat", "data_20190724_193122.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int), typemax(Int)]
const specs = [1:(length(amp) * length(det)), 1:(length(amp) * length(det))]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)

data_na = [datas_na[1]; datas_na[2]]

const prefix = joinpath(@__DIR__, "imgs", "data_20190724_na_img")

function show_ratio_image(xs, ys, data, idx; kws...)
    params, _ratios, uncs = NaCsData.get_values(data)
    @assert length(xs) * length(ys) == length(params)
    perm = sortperm(params)
    ratios = reshape(_ratios[perm, idx], length(xs), length(ys))
    X = zeros(length(xs) + 1)
    for i in 2:length(xs)
        X[i] = (xs[i - 1] + xs[i]) / 2
    end
    X[1] = 2 * xs[1] - X[2]
    X[end] = 2 * xs[end] - X[end - 1]
    Y = zeros(length(ys) + 1)
    for i in 2:length(ys)
        Y[i] = (ys[i - 1] + ys[i]) / 2
    end
    Y[1] = 2 * ys[1] - Y[2]
    Y[end] = 2 * ys[end] - Y[end - 1]
    pcolormesh(X, Y, ratios'; shading="flat", kws...)
end

show_loading_image(xs, ys, data; kws...) = show_ratio_image(xs, ys, data, 1; kws...)
show_survival_image(xs, ys, data; kws...) = show_ratio_image(xs, ys, data, 2; kws...)

figure()
show_survival_image(det, amp, data_na)
title("Na imaging survival")
xlabel("Detuning (MHz)")
ylabel("Amplitude")
colorbar()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
