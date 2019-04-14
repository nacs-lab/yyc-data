#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsPlot
using PyPlot
using LsqFit

function find_crossings(x, y)
    idxs = Float64[]
    for i in 1:(length(y) - 1)
        v1 = y[i]
        v2 = y[i + 1]
        if v1 <= 0 && v2 > 0
            push!(idxs, -v2 / (v1 - v2) * x[i] + v1 / (v1 - v2) * x[i + 1])
        end
    end
    return idxs
end

model_sin(x, p) = p[1] .* sin.(2 * π * p[2] .* x .+ p[3]) .+ p[4]
model_lin(x, p) = p[1] .* x .+ p[2]

function fit_sin(x, y)
    crossings = find_crossings(x, y)
    fit_zeros = curve_fit(model_lin, 0:(length(crossings) - 1), crossings,
                          [(crossings[end] - crossings[1]) / (length(crossings) - 1),
                           crossings[1]])
    maxy = maximum(y)
    miny = minimum(y)
    amp_est = (maxy - miny) / 2
    off_est = (maxy + miny) / 2
    f_est = 1 / fit_zeros.param[1]
    phase_est = ((-fit_zeros.param[2] / fit_zeros.param[1]) % 1) * 2π
    param_est = [amp_est, f_est, phase_est, off_est]
    fit = curve_fit(model_sin, x, y, param_est)
    param = fit.param
    unc = estimate_errors(fit)

    lo = minimum(x)
    hi = maximum(x)
    span = hi - lo
    mid = (hi + lo) / 2
    plot_lo = mid - span * 1.1 / 2
    plot_hi = mid + span * 1.1 / 2
    plotx = linspace(plot_lo, plot_hi, max(10000, length(x) * 10))
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model_sin.(plotx, (param,)))
end
