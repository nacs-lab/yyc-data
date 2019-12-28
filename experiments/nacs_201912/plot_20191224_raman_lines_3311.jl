#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

p = [8, 7.5, 4]
r = [367.405, 367.435, 367.601]
uncs = [0.02, 0.02, 0.039]

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

model_lin(x, p) = p[1] .+ p[2] .* x

fit = fit_data(model_lin, p, r, uncs, [367.8, -0.05])

const prefix = joinpath(@__DIR__, "imgs", "data_20191224_raman_lines")

figure()
errorbar(p, r, uncs, fmt="C0.")
plot(fit.plotx, fit.ploty, "C0")
text(4.3, 367.6, "\$f_{Raman}=f_{Raman0}+a\\cdot P\$")
text(4, 367.4, ("\$f_{Raman0}=$(fit.uncs[1]) \\mathrm{MHz}\$\n" *
                "\$a=$(fit.uncs[2]) \\mathrm{MHz\\cdot mW^{-1}}\$"),
     fontsize="small")
grid()
title("Raman light shift (306540.7 GHz)")
ylim([367.38, 367.66])
xlabel("Power (mW)")
ylabel("Resonance (MHz)")
NaCsPlot.maybe_save("$(prefix)_stark_power")

NaCsPlot.maybe_show()
