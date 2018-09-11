#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180908_090538.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((-1, 2,), (4,)))
data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, -2), (3,)))
data_nacs_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2,), (4,)))
data_nacs_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2,), (3,)))

const cspower = [5.0, 6.5, 8.3, 10.8, 13.9, 18.0, 23.2, 30.0]
const napower = [2.0, 2.8, 3.9, 5.4, 7.5, 10.4, 14.4, 20.0]

const spec_a = OrderedDict(
    :x0=>1:(length(cspower) * length(napower)),
    :x1=>1:(length(cspower) * length(napower)),
    :y0=>1:(length(cspower) * length(napower)),
    :y1=>1:(length(cspower) * length(napower)),
    :z0=>1:(length(cspower) * length(napower)),
    :z1=>1:(length(cspower) * length(napower)),
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_na_a = NaCsData.split_data(data_na_a, spec_a)
const split_nacs_cs_a = NaCsData.split_data(data_nacs_cs_a, spec_a)
const split_nacs_na_a = NaCsData.split_data(data_nacs_na_a, spec_a)

data_cs_0 = [split_cs_a[:x0]; split_cs_a[:y0]; split_cs_a[:z0]]
data_cs_x = split_cs_a[:x1]
data_cs_y = split_cs_a[:y1]
data_cs_z = split_cs_a[:z1]
data_na_0 = [split_na_a[:x0]; split_na_a[:y0]; split_na_a[:z0]]
data_na_x = split_na_a[:x1]
data_na_y = split_na_a[:y1]
data_na_z = split_na_a[:z1]
data_nacs_cs_0 = [split_nacs_cs_a[:x0]; split_nacs_cs_a[:y0]; split_nacs_cs_a[:z0]]
data_nacs_na_0 = [split_nacs_na_a[:x0]; split_nacs_na_a[:y0]; split_nacs_na_a[:z0]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180908_090538_merge_temp")

function show_survival_image(xs, ys, data; kws...)
    params, _ratios, uncs = NaCsData.get_values(data)
    @assert length(xs) * length(ys) == length(params)
    perm = sortperm(params)
    ratios = reshape(_ratios[perm, 2], length(xs), length(ys))
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
show_survival_image(cspower, napower, data_cs_0)
title("Cs / Cs")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_cs0")

figure()
show_survival_image(cspower, napower, data_na_0)
title("Na / Na")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_na0")

figure()
show_survival_image(cspower, napower, data_nacs_cs_0)
title("Cs / Na+Cs")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_nacs_cs0")

figure()
show_survival_image(cspower, napower, data_nacs_na_0)
title("Na / Na+Cs")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_nacs_na0")

figure()
show_survival_image(cspower, napower, data_cs_x)
title("Cs X cooling")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
show_survival_image(cspower, napower, data_cs_y)
title("Cs Y cooling")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_cs_ry")

figure()
show_survival_image(cspower, napower, data_cs_z)
title("Cs Z (axial) cooling")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
show_survival_image(cspower, napower, data_na_x)
title("Na X cooling")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
show_survival_image(cspower, napower, data_na_y)
title("Na Y cooling")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
show_survival_image(cspower, napower, data_na_z)
title("Na Z (axial) cooling")
xlabel("Cs power")
ylabel("Na power")
colorbar()
NaCsPlot.maybe_save("$(prefix)_na_az")

NaCsPlot.maybe_show()
