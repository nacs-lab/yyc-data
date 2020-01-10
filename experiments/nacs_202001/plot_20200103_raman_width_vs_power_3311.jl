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

powers = [2, 4, 8]
widths = [5.3, 10.4, 22.6]
uncs = [1.4, 3.7, 7.0]

model_pwr1(x, p) = p[1] .* x
model_pwr2(x, p) = p[1] .* x.^2
model_pwr(x, p) = p[1] .* x.^p[2]

fit_pwr1 = fit_data(model_pwr1, powers, widths, uncs, [2.6], plot_lo=0.0)
fit_pwr2 = fit_data(model_pwr2, powers, widths, uncs, [1.0], plot_lo=0.0)
fit_pwr = fit_data(model_pwr, powers, widths, uncs, [2.6, 1.0], plot_lo=0.0)

# @show fit_pwr1.uncs
# @show fit_pwr2.uncs
# @show fit_pwr.uncs

const prefix = joinpath(@__DIR__, "imgs", "data_20200103_raman_width_vs_power_3311")

figure()
errorbar(powers, widths, uncs, fmt="C0.")
plot(fit_pwr1.plotx, fit_pwr1.ploty, "C0", label="Linear")
plot(fit_pwr2.plotx, fit_pwr2.ploty, "C1", label="Quadratic")
plot(fit_pwr.plotx, fit_pwr.ploty, "C2", label="Power Law")
text(4.7, 5, "\$\\Gamma=a \\cdot P\$", color="C0")
text(3.3, 1.3, "\$a=$(fit_pwr1.uncs[1]) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$", color="C0")
text(0.1, 15, "\$\\Gamma=a \\cdot P^n\$", color="C2")
text(0.15, 8.6, ("\$a=$(fit_pwr.uncs[1])\$\n" *
                 "\$n=$(fit_pwr.uncs[2])\$"), color="C2")
legend(fontsize="small", loc="upper left")
grid()
title("Raman Linewidth (306480.7 GHz)")
xlim([0, 8.6])
ylim([0, 30])
xlabel("Tweezer Power (mW)")
ylabel("Raman Linewidth (kHz)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
