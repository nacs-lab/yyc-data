#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_1040 = joinpath(@__DIR__, "data", "sfg_wp_1040_20191029.csv")
const data_1040 = readdlm(iname_1040, ',', Float64, skipstart=1)
const iname_1550 = joinpath(@__DIR__, "data", "sfg_wp_1550_20191029.csv")
const data_1550 = readdlm(iname_1550, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "sfg_wp_20191029")

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

model(x, p) = p[1] .+ p[2] .* sind.((x .- p[3]) .* 2).^2

fit_1040 = fit_data(model, data_1040[:, 1], data_1040[:, 2],
                    [0.0, maximum(data_1040[:, 2]), data_1040[1, 1]])
figure()
plot(data_1040[:, 1], data_1040[:, 2], "C2o")
plot(fit_1040.plotx, fit_1040.ploty, "C0")
title("1040 Power calibration")
grid()
ylim([0, ylim()[2]])
text(-2, 5, "\$a+b\\cdot\\sin^2(2(\\theta - \\theta_0))\$", fontsize=18, color="C0")
text(17, 0.25, "\$a=$(fit_1040.uncs[1] * 1000)\$ mW\n\$b=$(fit_1040.uncs[2])\$ W\n\$\\theta_0=$(fit_1040.uncs[3])^\\circ\$", fontsize=20, color="C0")
xlabel("HWP Angle (deg)")
ylabel("Power (W)")
NaCsPlot.maybe_save("$(prefix)_1040")

fit_1550 = fit_data(model, data_1550[:, 1], data_1550[:, 2],
                    [0.0, maximum(data_1550[:, 2]), data_1550[1, 1]])
figure()
plot(data_1550[:, 1], data_1550[:, 2], "C2o")
plot(fit_1550.plotx, fit_1550.ploty, "C0")
title("1550 Power calibration")
grid()
ylim([0, ylim()[2]])
text(-2, 5, "\$a+b\\cdot\\sin^2(2(\\theta - \\theta_0))\$", fontsize=18, color="C0")
text(23, 0.25, "\$a=$(fit_1550.uncs[1] * 1000)\$ mW\n\$b=$(fit_1550.uncs[2])\$ W\n\$\\theta_0=$(fit_1550.uncs[3])^\\circ\$", fontsize=20, color="C0")
xlabel("HWP Angle (deg)")
ylabel("Power (W)")
NaCsPlot.maybe_save("$(prefix)_1550s")

NaCsPlot.maybe_show()
