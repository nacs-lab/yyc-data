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

const tweezer = [668.35, 565, 560, 555]
const freq = [1623, 616.54, 595.73, 577.41]
const freq_unc = [38, 0.38, 0.86, 0.25]

const prefix = joinpath(@__DIR__, "imgs", "data_20200402_ramann_vs_tweezer_3322")

function tweezer_model(x, p)
    return p[1] .- p[2] ./ (x .- p[3])
end

function gen_tweezer_model_res(f)
    return function(x, p)
        return p[1] .- p[2] ./ (x .- f)
    end
end

fit_freq = fit_data(tweezer_model, tweezer, freq, freq_unc, [250.0, 60000, 705],
                    plot_lo=300, plot_hi=670)
fit_freq2 = fit_data(gen_tweezer_model_res(705),
                     tweezer[2:end], freq[2:end], freq_unc[2:end], [250.0, 60000],
                     plot_lo=300, plot_hi=670)
@show fit_freq.uncs
@show fit_freq2.uncs

figure()
errorbar(tweezer, freq, freq_unc, fmt="C0s")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="With 668.35 GHz point")
plot(fit_freq2.plotx, fit_freq2.ploty, "C1", label="\$f_{res}=705\\ GHz\$")
legend(fontsize="small")
grid()
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("Raman Resonance (770XXX kHz)")
title("Raman Resonance")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
