#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200206_155743.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [-3.2 .+ linspace(-125, 125, 26)]
const specs_cs = [linspace(-80, 120, 26)]
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

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200206_155743_coprop_det")

function rabiLine(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end
function gen_rabi(t)
    return (x, p) -> p[1] .* rabiLine.(2π .* (x .- p[2]), t, p[3])
end
fit_na = fit_survival(gen_rabi(18.15e-3), datas_na[1], [0.9, 0, 100])
fit_cs = fit_survival(gen_rabi(16.88e-3), datas_cs[1], [0.9, 0, 100])

figure()
NaCsPlot.plot_survival_data(datas_na[1], fmt="C0.", label="Na")
plot(fit_na.plotx, fit_na.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_cs[1], fmt="C1.", label="Cs")
plot(fit_cs.plotx, fit_cs.ploty, "C1-")
title("Coprop Raman")
ylim([0, 1])
xlim([-152, 152])
text(-146, 0.78, ("\$\\Omega_{Na}=2\\pi\\!\\cdot\\!$(fit_na.uncs[3] / 2π)\$kHz\n" *
                 "\$f_{Na}=\\!$(fit_na.uncs[2])\$kHz\n" *
                 "\$p_{Na}=$(fit_na.uncs[1])\$"),
     color="C0", fontsize="x-small")
text(26, 0.78, ("\$\\Omega_{Cs}=2\\pi\\!\\cdot\\!$(fit_cs.uncs[3] / 2π)\$kHz\n" *
                "\$f_{Cs}=\\!$(fit_cs.uncs[2])\$kHz\n" *
                "\$p_{Cs}=$(fit_cs.uncs[1])\$"),
     color="C1", fontsize="x-small")
legend(fontsize="small", loc="center left")
grid()
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
