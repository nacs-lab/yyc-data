#!/usr/bin/julia -f

module Fitting

using LsqFit
using LinearAlgebra

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: linspace
import ..get_values, ..SortedData

export fit_data

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
        @static if VERSION >= v"1.0"
            fit = curve_fit(model, params, ratios, uncs.^-2, p0)
        else
            fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
        end
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = sqrt.(diag(estimate_covar(fit)))
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=plotx === false ? false : model.(plotx, (fit.param,)))
end

function fit_survival(model, data::SortedData, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end
end
