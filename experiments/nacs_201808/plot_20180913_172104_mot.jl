#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180913_172104.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))

const amp = linspace(0.2, 0.8, 16)
const det = -5:-1.0:-20

const prefix = joinpath(@__DIR__, "imgs", "data_20180913_172104_mot")

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
show_loading_image(amp, det, data_na_a)
title("Na MOT loading")
xlabel("Amplitude")
ylabel("Detuning (MHz)")
colorbar()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
