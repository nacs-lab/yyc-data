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

powers = [0.125, 0.25, 0.5, 1, 2, 4, 8]
freqs = [297.8415, 297.93787, 298.13086, 298.50836, 299.2481, 300.6769, 303.3651]
freq_uncs = [0.0011, 0.00067, 0.00068, 0.00085, 0.0012, 0.0038, 0.0096]
widths = [8.3, 6.8, 6.7, 14.6, 23.9, 79, 215]
width_uncs = [5.8, 3.7, 2.4, 3.5, 4.9, 15, 38]
# powers = [0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
# freqs = [297.8415, 297.93787, 298.13086, 298.50836, 299.2481, 300.6769, 303.3651, 308.217]
# freq_uncs = [0.0011, 0.00067, 0.00068, 0.00085, 0.0012, 0.0038, 0.0096, 0.23]
# widths = [8.3, 6.8, 6.7, 14.6, 23.9, 79, 215, 510]
# width_uncs = [5.8, 3.7, 2.4, 3.5, 4.9, 15, 38, 180]

model_pwr1(x, p) = p[1] .* (x)
model_pwr2(x, p) = p[1] .* (x).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit_freq = fit_data(model_lin_off, powers, freqs, freq_uncs, [297.8, 0.8],
                    plot_lo=0.0, plot_hi=9)
fit_freq2 = fit_data(model_sqr_off, powers, freqs, freq_uncs, [fit_freq.param; 0],
                     plot_lo=0.0, plot_hi=9)
fit_width = fit_data(model_pwr1, powers, widths, width_uncs, [15.2],
                     plot_lo=0.0, plot_hi=9)
fit_width2 = fit_data(model_pwr2, powers, widths, width_uncs, [10.0],
                      plot_lo=0.0, plot_hi=9)
# fit_freq = fit_data(model_lin_off, powers, freqs, freq_uncs, [297.8, 0.8],
#                     plot_lo=0.0, plot_hi=18)
# fit_freq2 = fit_data(model_sqr_off, powers, freqs, freq_uncs, [fit_freq.param; 0],
#                      plot_lo=0.0, plot_hi=18)
# fit_width = fit_data(model_pwr1, powers, widths, width_uncs, [15.2, 0],
#                      plot_lo=0.0, plot_hi=18)
# fit_width2 = fit_data(model_pwr2, powers, widths, width_uncs, [10.0, 0],
#                       plot_lo=0.0, plot_hi=18)

const prefix = joinpath(@__DIR__, "imgs", "data_20200209_raman_vs_power")

figure()
errorbar(powers, widths, width_uncs, fmt="C0.")
plot(fit_width.plotx, fit_width.ploty, "C0")
plot(fit_width2.plotx, fit_width2.ploty, "C2")
# text(0.2, 48, "\$\\Gamma=a \\cdot P\$", color="C0")
# text(0.35, 41, "\$a=$(fit_width.uncs[1]) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$", color="C0")
grid()
title("Raman Linewidth (306491 GHz)")
xlim([0, 8.5])
ylim([0, 250])
# xlim([0, 17])
# ylim([0, 700])
xlabel("Tweezer Power (mW)")
ylabel("Raman Linewidth (kHz)")
NaCsPlot.maybe_save("$(prefix)_w")

figure()
errorbar(powers, widths ./ powers .* 5.96, width_uncs ./ powers .* 5.96, fmt="C0o-")
grid()
title("Linewidth (306491 GHz)")
xlim([0, 8.5])
ylim([0, 190])
xlabel("Tweezer Power (mW)")
ylabel("Effective Excited Linewidth (MHz)")
NaCsPlot.maybe_save("$(prefix)_ew")

figure()
errorbar(powers, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="Linear")
plot(fit_freq2.plotx, fit_freq2.ploty, "C2", label="Quadratic")
text(0.1, 302.9, "\$f=f_0+a \\cdot P\$", color="C0", fontsize="small")
text(0.13, 301.6, ("\$f_0=$(fit_freq.uncs[1]) \\mathrm{MHz}\$\n" *
                   "\$a=$(fit_freq.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$"),
     color="C0", fontsize="small")
text(4.2, 299.6, "\$f=f_0+a \\cdot P+b \\cdot P^2\$", color="C2", fontsize="small")
text(3.3, 297.6, ("\$f_0=$(fit_freq2.uncs[1]) \\mathrm{MHz}\$\n" *
                  "\$a=$(fit_freq2.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$\n" *
                  "\$b=$(fit_freq2.uncs[3] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-2}\$"),
     color="C2", fontsize="small")
legend(fontsize="small", loc="upper left", ncol=2)
grid()
title("Raman Resonance (306491 GHz)")
xlim([0, 8.5])
# xlim([0, 17])
xlabel("Tweezer Power (mW)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_f")

figure()
errorbar(powers, (freqs .- model_lin_off(powers, fit_freq.param)) .* 1000,
         freq_uncs .* 1000, fmt="C0.-", label="Linear")
errorbar(powers, (freqs .- model_sqr_off(powers, fit_freq2.param)) .* 1000,
         freq_uncs .* 1000, fmt="C2.-", label="Quadratic")
legend(fontsize="small")
grid()
title("Frequency Fit Residue")
xlim([0, 8.5])
# xlim([0, 17])
ylim([-30, 30])
xlabel("Tweezer Power (mW)")
ylabel("Residue (kHz)")
NaCsPlot.maybe_save("$(prefix)_fres")

NaCsPlot.maybe_show()
