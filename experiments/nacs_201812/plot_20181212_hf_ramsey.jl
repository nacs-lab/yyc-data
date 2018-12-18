#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181212_162800.mat",
                "data_20181213_144356.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs_na = [0:2.0:200,
                  0:0.4:40]
const specs_cs = [0:2.0:200,
                  0:2.0:200]
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

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na = datas_na[1]
data_cs_raman = datas_cs[1]
data_cs_uwave = datas_cs[2]

model_cosexp(x, p) = p[1] .+ p[2] .* cos.(2π .* x .* p[3] .+ p[4]) .* exp.(.-x ./ p[5])

fit_na = fit_survival(model_cosexp, data_na, [0.5, 0.4, 0.01, 2, 150])
@show fit_na.uncs
fit_cs_raman = fit_survival(model_cosexp, data_cs_raman, [0.5, 0.4, 0.025, 1, 100])
@show fit_cs_raman.uncs
fit_cs_uwave = fit_survival(model_cosexp, data_cs_uwave, [0.5, 0.4, 0.019, 1, 100])
@show fit_cs_uwave.uncs

const prefix = joinpath(@__DIR__, "imgs", "data_20181212_hf_ramsey")

figure()
NaCsPlot.plot_survival_data(data_na, fmt="C0.", label="Na")
plot(fit_na.plotx, fit_na.ploty, "C0-")
NaCsPlot.plot_survival_data(data_cs_raman, fmt="C1.", label="Cs Raman")
plot(fit_cs_raman.plotx, fit_cs_raman.ploty, "C1-")
NaCsPlot.plot_survival_data(data_cs_uwave, fmt="C2.", label="Cs \$\\mu\$Wave")
plot(fit_cs_uwave.plotx, fit_cs_uwave.ploty, "C2-")
text(25, 0.05, "\$\\tau_{Na}=$(fit_na.uncs[5])\\mu s\$", color="C0", fontsize="small")
text(0, 0.9, "\$\\tau_{Cs}^{Raman}=$(fit_cs_raman.uncs[5])\\mu s\$",
     color="C1", fontsize="small")
text(115, 0.1, "\$\\tau_{Cs}^{\\mu Wave}=$(fit_cs_uwave.uncs[5])\\mu s\$",
     color="C2", fontsize="small")
grid()
legend(fontsize="small")
ylim([0, 1])
xlim([0, 220])
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
