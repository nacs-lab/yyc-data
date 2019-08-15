#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190814_225957.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 ]
const specs = [([0.0, 1, 2, 5, 10, 20, 50, 80, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 80, 100])]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function fit_survival(model, data, p0; plotx=nothing, plot_lo=nothing, plot_hi=nothing,
                      use_unc=true, plot_scale=1.1)
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
        fit = curve_fit(model, params, ratios[:, 2], uncs[:, 2].^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

model_exp1(x, p) = p[1] .* exp.(x ./ -p[2])
function gen_model2(p0)
    model(x, p) = p0[1] .* exp.(x ./ -p0[2]) .* (p[1] .* exp.(x ./ -p[2]) .+ p[3])
end

data_nopa = datas_nacs[1][1]
data_pa = datas_nacs[1][2]

fit_nopa = fit_survival(model_exp1, data_nopa, [0.9, 500])
fit_pa = fit_survival(gen_model2(fit_nopa.param), data_pa, [0.6, 10, 0.4])

const prefix = joinpath(@__DIR__, "imgs", "data_20190814_225957_pa_v14_rate")

figure()
NaCsPlot.plot_survival_data(data_nopa, fmt="C0.", label="No PA")
plot(fit_nopa.plotx, fit_nopa.ploty, "C0-")
NaCsPlot.plot_survival_data(data_pa, fmt="C1.", label="PA")
plot(fit_pa.plotx, fit_pa.ploty, "C1-")
legend()
grid()
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_xticks]([1, 10, 100])
ax[:set_xticklabels](["1", "10", "100"])
ylim([0.3, 0.9])
xlim([0.5, 110])
text(3, 0.8, "\$\\tau_{NoPA}=$(fit_nopa.uncs[2] / 1000) s\$", color="C0")
text(1.2, 0.6, "\$\\tau_{PA}=$(fit_pa.uncs[2]) ms\$", color="C1")
title("v=14 PA")
xlabel("PA time (ms)")
ylabel("Two body survival")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
