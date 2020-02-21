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

tweezer = [460.7, 480.7, 490.7, 522, 530, 550.7]
widths = [17.6, 22.6, 52, 430, 118, 20.7]
width_uncs = [5, 7, 10, 110, 33, 4.5]
# tweezer = [460.7, 480.7, 490.7, 530, 550.7]
# widths = [17.6, 22.6, 52, 118, 20.7]
# width_uncs = [5, 7, 10, 33, 4.5]

const prefix = joinpath(@__DIR__, "imgs", "data_20200220_raman_vs_freq_3311")

figure()
errorbar(tweezer, widths, width_uncs, fmt="C0s-")
grid()
title("Raman Linewidth (8 mW)")
xlabel("Tweezer Frequency (306XXX GHz)")
ylabel("Raman Linewidth (kHz)")
ylim([0, 150])
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
