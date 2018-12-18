#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181217_220313.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_cs = [OrderedDict(
    :f4=>0:8.0:200,
    :f14=>0:8.0:200,
    :f24=>0:8.0:200
)]
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

const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_f4 = datas_cs[1][:f4]
data_f14 = datas_cs[1][:f14]
data_f24 = datas_cs[1][:f24]

const prefix = joinpath(@__DIR__, "imgs", "data_20181217_220313_uwave_ramsey_freq")

model_cosexp(x, p) = p[1] .+ p[2] .* cos.(2π .* x .* p[3] .+ p[4]) .* exp.(.-x ./ p[5])

fit_f4 = fit_survival(model_cosexp, data_f4, [0.5, 0.4, 0.02, 0, 100])
@show fit_f4.uncs
fit_f14 = fit_survival(model_cosexp, data_f14, [0.5, 0.4, 0.008, 0, 100])
@show fit_f14.uncs
fit_f24 = fit_survival(model_cosexp, data_f24, [0.5, 0.4, 0.005, 0, 100])
@show fit_f24.uncs

figure()
NaCsPlot.plot_survival_data(data_f4, fmt="C0.", label="4kHz")
plot(fit_f4.plotx, fit_f4.ploty, "C0-")
NaCsPlot.plot_survival_data(data_f14, fmt="C1.", label="14kHz")
plot(fit_f14.plotx, fit_f14.ploty, "C1-")
NaCsPlot.plot_survival_data(data_f24, fmt="C2.", label="24kHz")
plot(fit_f24.plotx, fit_f24.ploty, "C2-")
text(40, 0.9, "\$\\nu_{4}=$(fit_f4.uncs[3] * 1000)kHz\$", color="C0")
text(50, 0.02, "\$\\nu_{14}=$(fit_f14.uncs[3] * 1000)kHz\$", color="C1")
text(100, 0.12, "\$\\nu_{24}=$(fit_f24.uncs[3] * 1000)kHz\$", color="C2")
grid()
legend(fontsize="small")
ylim([0, 1])
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
