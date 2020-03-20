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

Bs = [100, 95, 80, 60, 40]
freqs = [770.59514, 770.5767, 770.521, 770.4479, 770.3713]
freq_uncs = [0.00054, 0.0018, 0.004, 0.0027, 0.0031]

model_pwr1(x, p) = p[1] .* (x)
model_pwr2(x, p) = p[1] .* (x).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit_freq = fit_data(model_lin_off, Bs, freqs, freq_uncs, [770.2, 0.003])
# fit_freq2 = fit_data(model_sqr_off, Bs, freqs, freq_uncs, [770.2, 25, 0])

const prefix = joinpath(@__DIR__, "imgs", "data_20200317_raman_3322_vs_B")

figure()
errorbar(Bs, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="Linear")
# plot(fit_freq2.plotx, fit_freq2.ploty, "C2", label="Quadratic")
text(62, 770.41, "\$f=f_0+a \\cdot B\$", color="C0", fontsize="small")
text(55, 770.355, ("\$f_0=$(fit_freq.uncs[1]) \\mathrm{MHz}\$\n" *
                   "\$a=$(fit_freq.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{\\%}^{-1}\$"),
     color="C0", fontsize="small")
# text(9, 770.38, "\$f=f_0+a \\cdot B+b \\cdot B^2\$", color="C2", fontsize="small")
# text(7, 770.27, ("\$f_0=$(fit_freq2.uncs[1]) \\mathrm{MHz}\$\n" *
#                  "\$a=$(fit_freq2.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{\\%}^{-1}\$\n" *
#                  "\$b=$(fit_freq2.uncs[3] * 1000) \\mathrm{kHz}\\cdot\\mathrm{\\%}^{-2}\$"),
#      color="C2", fontsize="small")
# legend(fontsize="small", loc="upper left", ncol=2)
grid()
title("Raman Resonance 288560 GHz 15 mW")
xlabel("B field (%)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(Bs, (freqs .- model_lin_off(Bs, fit_freq.param)) .* 1000,
         freq_uncs .* 1000, fmt="C0.-", label="Linear")
# errorbar(Bs, (freqs .- model_sqr_off(Bs, fit_freq2.param)) .* 1000,
#          freq_uncs .* 1000, fmt="C2.-", label="Quadratic")
# legend(fontsize="small")
grid()
title("Frequency Fit Residue")
xlabel("B field (%)")
ylabel("Residue (kHz)")
NaCsPlot.maybe_save("$(prefix)_res")

NaCsPlot.maybe_show()
