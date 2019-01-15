#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190113_191153.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [(0:1.0:20,
                   0:1.0:20)]
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

model_cos_off(x, p) = p[1] .+ (p[2] / 2) .* (cos.(x ./ p[3] .* π) .- 1)

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)

data_na_11 = datas_na[1][1]
data_na_10 = datas_na[1][2]

fit_na_11 = fit_survival(model_cos_off, data_na_11, [0.9, 0.8, 13])
fit_na_10 = fit_survival(model_cos_off, data_na_10, [0.9, 0.1, 11])

const prefix = joinpath(@__DIR__, "imgs", "data_20190113_191153_depump_coprop_time")

figure()
NaCsPlot.plot_survival_data(data_na_11, fmt="C0.", label="1, 1")
plot(fit_na_11.plotx, fit_na_11.ploty, "C0")
NaCsPlot.plot_survival_data(data_na_10, fmt="C1.", label="1, 0")
plot(fit_na_10.plotx, fit_na_10.ploty, "C1")
legend(fontsize="small", handlelength=1, handleheight=0.4,
       handletextpad=0.3, loc="lower left")
text(7, 0.5, "\$\\tau_{\\pi\\ 1, 1}=$(fit_na_11.uncs[3])\\mu s\$", color="C0")
text(9, 0.05, "\$p_{1, 1}=$(fit_na_11.uncs[2])\$", color="C0")
text(7, 0.9, "\$\\tau_{\\pi\\ 1, 0}=$(fit_na_10.uncs[3])\\mu s\$", color="C1")
text(8, 0.7, "\$p_{1, 0}=$(fit_na_10.uncs[2])\$", color="C1")
grid()
ylim([0, 1])
title("Na depump state distribution")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

NaCsPlot.maybe_show()
