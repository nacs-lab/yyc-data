#!/usr/bin/julia -f

module Fitting

using LsqFit
using LinearAlgebra

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: linspace
import ..get_values, ..SortedData

export fit_data

function curve_fit(model, params::AbstractArray{T}, args...; kws...) where T<:Number
    return LsqFit.curve_fit(model, params, args...; kws...)
end

function curve_fit(origin_model, origin_xs, args...; kws...)
    xs = 1:length(origin_xs)
    function model(x, p)
        x == xs && return origin_model(origin_xs, p)
        return origin_model(origin_xs[x], p)
    end
    return LsqFit.curve_fit(model, xs, args...; kws...)
end

function get_plot_range(xs; plotx=nothing, plot_lo=nothing,
                        plot_hi=nothing, plot_scale=1.1, plot_npts=10000)
    if plotx !== nothing
        return plotx
    end
    lo = minimum(xs)
    hi = maximum(xs)
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
    plotx = linspace(plot_lo, plot_hi, plot_npts)
end

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, xs, ratios, uncs, p0; plotx=nothing, plot_lo=nothing,
                  plot_hi=nothing, plot_scale=1.1, plot_npts=10000, kws...)
    use_unc = uncs !== nothing
    if use_unc
        fit = curve_fit(model, xs, ratios, uncs.^-2, p0; kws...)
    else
        fit = curve_fit(model, xs, ratios, p0; kws...)
    end
    plotx = get_plot_range(xs; plotx=plotx, plot_lo=plot_lo,
                           plot_hi=plot_hi, plot_scale=plot_scale, plot_npts=plot_npts)
    param = fit.param
    unc = sqrt.(diag(estimate_covar(fit)))
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=plotx === false ? false : model(plotx, param))
end

function fit_survival(model, data::SortedData, p0; use_unc=true, kws...)
    if use_unc
        xs, ratios, uncs = get_values(data)
        return fit_data(model, xs, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        xs, ratios, uncs = get_values(data, 0.0)
        return fit_data(model, xs, ratios[:, 2], p0; kws...)
    end
end
end
