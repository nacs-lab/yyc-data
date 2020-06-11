#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const na_amp = linspace(0.3, 1, 8)
const na_det = -24:2.0:-12
const cs_amp = linspace(0.3, 1, 8)
const cs_det = -11:1.0:-5

const inames = ["data_20200610_211615.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [1:(length(na_amp) * length(na_det))]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

const data_na = datas_na[1]
const data_cs = datas_cs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20200610_211615_mot")

function show_loading_image(xs, ys, data; kws...)
    params, _ratios, uncs = NaCsData.get_values(data)
    @assert length(xs) * length(ys) == length(params)
    perm = sortperm(params)
    ratios = reshape(_ratios[perm, 1], length(xs), length(ys))
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

figure()
show_loading_image(na_amp, na_det, data_na)
title("Na CMOT loading")
xlabel("Amplitude")
ylabel("Detuning (MHz)")
colorbar()
NaCsPlot.maybe_save("$(prefix)_na")

figure()
show_loading_image(cs_amp, cs_det, data_cs)
title("Cs MOT loading")
xlabel("Amplitude")
ylabel("Detuning (MHz)")
colorbar()
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
