#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181010_183350.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const spec_a = linspace(-21, 21, 15) .+ 29

function select_times(selector, binsz)
    nparam = length(params_a)
    return [NaCsData.split_data(NaCsData.select_count(params_a, logicals_a, selector,
                                                      i->i in i0:(i0 + binsz - 1)),
                                spec_a) for i0 in 1:binsz:(nparam - binsz ÷ 2)]
end

const data_nacs_a = select_times(NaCsData.select_single((1, 2,), (3, 4,)), 1620 * 5)
const data_nacs_all = NaCsData.split_data(
    NaCsData.select_count(params_a, logicals_a,
                          NaCsData.select_single((1, 2,), (3, 4,))), spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20181010_183350_raman_rt")

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

function fit_survival(model, data, p0; plotx=nothing, use_unc=false, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        plotx = linspace(mid - span * plot_scale / 2, mid + span * plot_scale / 2, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2].^2, p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

raman_model(x, p) = p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ p[4]).^2)
fit = fit_survival(raman_model, data_nacs_all, [0.3, 0.3, 16, 10])

figure()
title("NaCs -> NaCs")
NaCsPlot.plot_survival_data(data_nacs_all, fmt="C0.")
plot(fit.plotx, fit.ploty, "C0-")
grid()
xlabel("Frequency (299300 + kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_all")

xs, img = get_time_img(data_nacs_a)
figure()
imshow(img, interpolation="gaussian", extent=[xs[1], xs[end], size(img, 1), 1];
       cmap="viridis_r", aspect="auto")
title("NaCs -> NaCs")
xlabel("Frequency (299300 + kHz)")
ylabel("Time index")
colorbar()
NaCsPlot.maybe_save("$(prefix)_time")

img2 = img .- raman_model(xs, fit.param)'
figure()
imshow(img2, interpolation="gaussian", extent=[xs[1], xs[end], size(img2, 1), 1];
       cmap="viridis_r", aspect="auto")
title("NaCs -> NaCs")
xlabel("Frequency (299300 + kHz)")
ylabel("Time index")
colorbar()
NaCsPlot.maybe_save("$(prefix)_diff")

NaCsPlot.maybe_show()
