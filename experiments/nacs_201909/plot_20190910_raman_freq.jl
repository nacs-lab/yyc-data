#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190910_231236.mat",
                "data_20190912_211929.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [732.2 .+ (-1:0.2:1),
               731.7 .+ linspace(-0.5, 0.5, 11)]
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

const prefix = joinpath(@__DIR__, "imgs", "data_20190910_raman_freq")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
fit1 = fit_survival(model_lorentzian, datas_nacs[1], [0.7, 0.3, 732, 0.5])
fit2 = fit_survival(model_lorentzian, datas_nacs[2], [0.7, 0.3, 731.6, 0.5])

figure()
NaCsPlot.plot_survival_data(datas_nacs[1], fmt="C0.", label="0910")
plot(fit1.plotx, fit1.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_nacs[2], fmt="C1.", label="0912")
plot(fit2.plotx, fit2.ploty, "C1-")
text(732.25, 0.485, "\$\\Gamma_{0910}=$(fit1.uncs[4])\$ kHz", color="C0", fontsize="small")
text(732.1, 0.46, "\$f_{0910}=$(fit1.uncs[3])\$ kHz", color="C0", fontsize="small")
text(731.5, 0.725, "\$\\Gamma_{0912}=$(fit2.uncs[4])\$ kHz", color="C1", fontsize="small")
text(731.46, 0.7, "\$f_{0912}=$(fit2.uncs[3])\$ kHz", color="C1", fontsize="small")
legend(fontsize="small", loc="lower left")
grid()
xlabel("297XXX kHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
