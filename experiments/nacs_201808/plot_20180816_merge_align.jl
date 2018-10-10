#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

names = ["data_20180816_233919.mat",
         "data_20180817_010250.mat",
         "data_20180817_022619.mat",
         "data_20180817_034950.mat",
         "data_20180817_051319.mat",
         "data_20180817_063649.mat",
         "data_20180817_080024.mat",
         "data_20180817_092412.mat",
         "data_20180817_104735.mat",
         "data_20180817_121112.mat",
         "data_20180817_133445.mat",
         "data_20180817_160117.mat",
         "data_20180817_172500.mat",
         "data_20180817_184838.mat",
         "data_20180817_201215.mat",
         "data_20180817_214528.mat"]

param_logicals = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", name))
                  for name in names]

data_nas = [NaCsData.select_count(p[1], p[2], NaCsData.select_single((1,), (3,)))
            for p in param_logicals]
data_css = [NaCsData.select_count(p[1], p[2], NaCsData.select_single((2,), (4,)))
            for p in param_logicals]

const spec = OrderedDict(
    :az=>2.9 .+ linspace(-0.6, 0.6, 13),
    :rx=>2.9 .+ linspace(-0.6, 0.6, 13),
    :ry=>2.9 .+ linspace(-0.6, 0.6, 13),
    :surv=>2.9 .+ linspace(-0.6, 0.6, 13),
    :align=>2.9 .+ linspace(-0.6, 0.6, 13),
)

split_nas = [NaCsData.split_data(d, spec) for d in data_nas]
split_css = [NaCsData.split_data(d, spec) for d in data_css]

const prefix = joinpath(@__DIR__, "imgs", "data_20180816_merge_align")
const data_prefix = joinpath(@__DIR__, "plot_data", "data_20180816_merge_align")

const ypos = (1:length(names)) .* 3 * 0.03

function show_survival_image(ys, datas; kws...)
    @assert length(ys) == length(datas)
    img = Matrix{Float64}(length(ys), size(datas[1], 1))
    local xs
    for i in 1:length(datas)
        params, ratios, uncs = NaCsData.get_values(datas[i])
        perm = sortperm(params)
        if !@isdefined(xs)
            xs = params[perm]
        end
        img[i, :] = ratios[perm, 2]
    end
    imshow(img, interpolation="nearest", extent=[xs[1], xs[end], ys[end], ys[1]]; kws...)
    data = Matrix{Float64}(length(ys) + 1, length(xs) + 1)
    data[1, 1] = 0
    data[2:end, 1] = ys
    data[1, 2:end] = xs
    data[2:end, 2:end] = img
    return data
end

figure()
data = show_survival_image(ypos, [split_css[i][:align] for i in 1:length(names)])
colorbar()
title("Cs alignment")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_cs_align.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_cs_align")

figure()
data = show_survival_image(ypos, [split_css[i][:az] for i in 1:length(names)], vmin=0)
colorbar()
title("Cs axial")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_cs_az.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
data = show_survival_image(ypos, [split_css[i][:rx] for i in 1:length(names)], vmin=0)
colorbar()
title("Cs radial X")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_cs_rx.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
data = show_survival_image(ypos, [split_css[i][:ry] for i in 1:length(names)], vmin=0)
colorbar()
title("Cs radial Y")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_cs_ry.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_cs_ry")

figure()
data = show_survival_image(ypos, [split_nas[i][:az] for i in 1:length(names)], vmin=0)
colorbar()
title("Na axial")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_na_az.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
data = show_survival_image(ypos, [split_nas[i][:rx] for i in 1:length(names)], vmin=0)
colorbar()
title("Na radial X")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_na_rx.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
data = show_survival_image(ypos, [split_nas[i][:ry] for i in 1:length(names)], vmin=0)
colorbar()
title("Na radial Y")
xlabel("dMerge (um)")
ylabel("Vertical position (um)")
writedlm("$(data_prefix)_na_ry.csv", data, ',')
NaCsPlot.maybe_save("$(prefix)_na_ry")

NaCsPlot.maybe_show()
