#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170918_215351.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const tlowers = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
const vlowers = [0.02, 0.04, 0.08, 0.1, 0.12]
const nt = length(tlowers)
const nv = length(vlowers)

const spec_a = OrderedDict(
    :base=>[0],
    :lower=>1:(nt * nv)
)

const split_a = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_a), spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20170918_215351_lower_2v")

@show split_a[:base]

params, ratios, uncs = NaCsData.get_values(split_a[:lower])
ratios = reshape(ratios[:, 2], (nt, nv))
uncs = reshape(uncs[:, 2], (nt, nv))

figure(figsize=[6, 8])
imshow(ratios)
title("Axial sideband")
xlabel("VLower (V)")
ylabel("TLower (ms)")
xticks(0:(nv - 1), vlowers)
yticks(0:(nt - 1), tlowers)
colorbar()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
