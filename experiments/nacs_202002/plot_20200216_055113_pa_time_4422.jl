#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200216_055113.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [([0, 1, 2, 5, 10, 25, 50, 100, 200] .* 1.0,
                [0, 1, 2, 5, 10, 25, 50, 100, 200] .* 2.0)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, params, ratios, uncs, p0;
                  plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
    use_unc = uncs !== nothing
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        if plot_lo === nothing
            plot_lo = mid - span * plot_scale / 2
            if plot_lo * lo <= 0
                plot_lo = 0
            end
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
            if plot_hi * hi <= 0
                plot_hi = 0
            end
        end
        plotx = linspace(plot_lo, plot_hi, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200216_055113_raman_time_4422")

function model_expoff(x, p)
    p[3] .+ p[1] .* exp.(.- x ./ p[2])
end
function model_exp(x, p)
    p[1] .* exp.(.- x ./ p[2])
end
fit1 = fit_survival(model_expoff, datas_nacs[1][1], [0.6, 150, 0.1])
fit2 = fit_survival(model_expoff, datas_nacs[1][2], [0.6, 600, 0.1])

# @show fit1.uncs
# @show fit2.uncs

figure()
NaCsPlot.plot_survival_data(datas_nacs[1][1], fmt="C0.", label="16 mW")
plot(fit1.plotx, fit1.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_nacs[1][2], fmt="C1.", label="8 mW")
plot(fit2.plotx, fit2.ploty, "C1-")
xlim([0, 410])
text(7.5, 0.065, "\$\\tau_{16}=$(fit1.uncs[2])\$ ms", color="C0")
text(100, 0.36, "\$\\tau_{8}=$(fit2.uncs[2])\$ ms", color="C1")
legend(fontsize="small", loc="upper right")
title("306492 GHz")
grid()
xlabel("Time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
