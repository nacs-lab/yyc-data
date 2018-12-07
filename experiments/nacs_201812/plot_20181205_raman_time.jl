#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181205_170215.mat",
                "data_20181205_193552.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [[0, 0.2, 0.4, 0.7, 1, 2, 4, 7],
               [0, 0.1, 0.2, 0.35, 0.5, 1, 2, 3.5]]

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
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
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

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2,), (3, 4,)), maxcnts, specs)
const data_r50 = datas_nacs[2]
const data_r200 = datas_nacs[1]

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

fit_r50 = fit_survival(model_exp_off, data_r50, [0.6, 0.8, 0.1], plot_lo=0)
fit_r200 = fit_survival(model_exp_off, data_r200, [0.6, 3.2, 0.1], plot_lo=0)

const prefix = joinpath(@__DIR__, "imgs", "data_20181205_raman_time")

figure()
NaCsPlot.plot_survival_data(data_r50, fmt="C0.", label="0.02")
plot(fit_r50.plotx, fit_r50.ploty, "C0-")
NaCsPlot.plot_survival_data(data_r200, fmt="C1.", label="0.005")
plot(fit_r200.plotx, fit_r200.ploty, "C1-")
grid()
text(3, 0.5, "\$a+b\\cdot e^{-\\dfrac{t}{\\tau}}\$", fontsize="large")
text(1.1, 0.21, "\$a=$(fit_r50.uncs[3])\$\n" *
     "\$b=$(fit_r50.uncs[1])\$\n" *
     "\$\\tau=$(fit_r50.uncs[2])\$ms", color="C0", fontsize="small")
text(4, 0.32, "\$a=$(fit_r200.uncs[3])\$\n" *
     "\$b=$(fit_r200.uncs[1])\$\n" *
     "\$\\tau=$(fit_r200.uncs[2])\$ms", color="C1", fontsize="small")
legend(fontsize="small", labelspacing=0.2, borderpad=0.2,
       handletextpad=0.2, columnspacing=0.3, borderaxespad=0.2)
xlabel("Raman time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
