#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190521_232054.mat",
                "data_20190522_073239.mat",
                "data_20190522_130331.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 ]
const specs = [[0.0, 2, 5, 10, 20, 50, 100],
               [0.0 0.2 0.5 1 2 4 6 8 10 20],
               [0.0 0.2 0.5 1 2 4 6 8 10 20]]
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

data_nacs_4 = datas_nacs[1]
data_nacs_3 = [datas_nacs[2]; datas_nacs[3]][[1; 3:6; 8:10]]

fit_nacs_4 = fit_survival(model_exp_off, data_nacs_4, [0.45, 20, 0.2])
fit_nacs_3 = fit_survival(model_exp_off, data_nacs_3, [0.45, 2, 0.2])

const prefix = joinpath(@__DIR__, "imgs", "data_20190521_lifetime_cmp")

figure()
NaCsPlot.plot_survival_data(data_nacs_4, fmt="C0.", label="4, 4")
plot(fit_nacs_4.plotx, fit_nacs_4.ploty, "C0-")
NaCsPlot.plot_survival_data(data_nacs_3, fmt="C1.", label="3, 3")
plot(fit_nacs_3.plotx, fit_nacs_3.ploty, "C1-")
text(0.15, 0.64, "\$\\tau_{4,4}=$(fit_nacs_4.uncs[2]) ms\$", color="C0")
text(0.3, 0.25, "\$\\tau_{3,3}=$(fit_nacs_3.uncs[2]) ms\$", color="C1")
legend()
grid()
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_xticks]([0.1, 1, 10, 100])
ax[:set_xticklabels](["0.1", "1", "10", "100"])
ylim([0.15, 0.7])
xlim([0.1, 100])
title("Photo Association")
xlabel("PA time (ms)")
ylabel("Na + Cs survival")
NaCsPlot.maybe_save("$(prefix)_log")

figure()
NaCsPlot.plot_survival_data(data_nacs_4, fmt="C0.", label="4, 4")
plot(fit_nacs_4.plotx, fit_nacs_4.ploty, "C0-")
NaCsPlot.plot_survival_data(data_nacs_3, fmt="C1.", label="3, 3")
plot(fit_nacs_3.plotx, fit_nacs_3.ploty, "C1-")
text(4, 0.60, "\$\\tau_{4,4}=$(fit_nacs_4.uncs[2]) ms\$", color="C0")
text(4, 0.55, "\$\\tau_{3,3}=$(fit_nacs_3.uncs[2]) ms\$", color="C1")
legend()
grid()
ylim([0.25, 0.7])
xlim([0, 30])
title("Photo Association")
xlabel("PA time (ms)")
ylabel("Na + Cs survival")
NaCsPlot.maybe_save("$(prefix)_lin")

NaCsPlot.maybe_show()
