#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190802_202451.mat",
                "data_20190803_002700.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 ]
const specs = [[1.0, 2, 3, 4, 5, 8, 10, 20, 50],
               ([1.0, 2, 3, 4, 5, 8, 10, 20, 50],
                [1.0, 2, 3, 4, 5, 8, 10, 20, 50])]
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

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

data_nopa = datas_nacs[2][1]
data_pa = [datas_nacs[1]; datas_nacs[2][2]]

fit_nopa = fit_survival(model_exp_off, data_nopa, [0.3, 3, 0])
fit_pa = fit_survival(model_exp_off, data_pa, [0.3, 10, 0])
τ_paonly = 1 / (1 / fit_pa.uncs[2] - 1 / fit_nopa.uncs[2])

const prefix = joinpath(@__DIR__, "imgs", "data_20190802_pa_rate_33_22")

figure()
NaCsPlot.plot_survival_data(data_nopa, fmt="C0.", label="No PA")
plot(fit_nopa.plotx, fit_nopa.ploty, "C0-")
NaCsPlot.plot_survival_data(data_pa, fmt="C1.", label="PA")
plot(fit_pa.plotx, fit_pa.ploty, "C1-")
legend()
grid()
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_xticks]([1, 10])
ax[:set_xticklabels](["1", "10"])
ylim([0, 0.35])
xlim([0.5, 60])
text(0.6, 0.28, "\$\\tau_{NoPA}=$(fit_nopa.uncs[2])\$", color="C0")
text(0.6, 0.01, "\$\\tau_{PA}=$(fit_pa.uncs[2])\$", color="C1")
text(1.2, 0.15, "\$\\tau_{\\Delta}=$(τ_paonly)\$", color="C3")
title("3, 3 + 2, 2 PA")
xlabel("PA time (ms)")
ylabel("Two body survival")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
