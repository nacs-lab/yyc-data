#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180421_002556.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if !logicals[1] || !logicals[2]
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3] && logicals[4]]
end

data_2body = NaCsData.select_count(params_a, logicals_a, selector)

const spec_a = (298.5 .+ linspace(-10, 10, 51),
                298.5 .+ linspace(-5, 5, 51),
                298.5 .+ linspace(-3, 3, 51),
                298.5 .+ linspace(-3, 3, 51),
                298.5 .+ linspace(-3, 3, 51))

const split_2body = NaCsData.split_data(data_2body, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180421_002556")

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

function eit_model(x, p)
    return p[1] .+ p[2] ./ (1 .+ (x .- p[3]).^2 ./ (p[4] ./ 2).^2)
end

const fit_2body_2 = fit_survival(eit_model, split_2body[2], [0.3, 0.45, 298.4, 1.5])

figure()
NaCsPlot.plot_survival_data(split_2body[2], fmt="C0o")
plot(fit_2body_2.plotx, fit_2body_2.ploty, "k")
grid()
title("Dark resonance spectrum")
xlabel("Two photon detuning (MHz)")
ylabel("Probability of both survive")
NaCsPlot.maybe_save("$(prefix)_eit-v0")

NaCsPlot.maybe_show()
