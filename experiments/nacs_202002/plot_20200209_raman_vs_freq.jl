#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

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

tweezer = [480, 485, 490, 491]
freqs = [299.5381, 300.4961, 304.595, 308.241]
freq_uncs = [0.0016, 0.0053, 0.034, 0.045]

model_det(x, p) = p[1] .- p[2] ./ (x .- p[3])
model_pwr1(x, p) = p[1] .* (x .+ p[2])
model_pwr2(x, p) = p[1] .* (x .+ p[2]).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit_freq = fit_data(model_det, tweezer, freqs, freq_uncs, [298, 20, 493.0])

const prefix = joinpath(@__DIR__, "imgs", "data_20200209_raman_vs_freq")

figure()
errorbar(tweezer, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0")
text(480.2, 309, "\$f_{Raman}=f_{Raman0}-\\dfrac{a}{f_{PA}-f_{PA0}}\$",
     fontsize="small", color="C0")
text(480.7, 304.2, ("\$f_{Raman0}=$(fit_freq.uncs[1])\\ \\mathrm{MHz}\$\n" *
                    "\$a=$(fit_freq.uncs[2])\\ \\mathrm{MHz\\cdot GHz}\$\n" *
                    "\$f_{PA0}=$(fit_freq.uncs[3] + 306000)\\ \\mathrm{GHz}\$"),
     fontsize="small", color="C0")
grid()
title("Raman Resonance (16 mW)")
xlabel("Tweezer Frequency (306XXX GHz)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
