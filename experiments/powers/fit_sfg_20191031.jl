#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci

const prefix = joinpath(@__DIR__, "imgs", "sfg_20191031")

P1550(θ) = 7.43 .* sind.(2 .* (θ - 2.44)).^2
P1040(θ) = 9.94 .* sind.(2 .* (θ + 1.4)).^2

const Input1550 = [0.150; P1550.([8, 12, 15, 18, 26, 36, 47])]
const Input1040 = [0.150; P1040.([4, 6, 10, 16, 24, 32, 44])]
const Output = [0.00252, 0.01030, 0.0440, 0.165,
                0.451, 1.70 / 0.479 * 0.451, 3.44 / 0.479 * 0.451, 4.70 / 0.479 * 0.451]

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

model(x, p) = p[1] .* x

fit_output = fit_data(model, Input1040[1:3] .* Input1550[1:3],
                      Output[1:3], [100.0], plot_hi=Input1040[end] .* Input1550[end] * 1.3)

figure()
plot(fit_output.plotx, fit_output.ploty, "C1-")
plot(Input1040 .* Input1550, Output, "C0o")
grid()
text(0.2, 0.003, "\$$(fit_output.uncs[1])\\mathrm{W}^{-1}\$", fontsize=18, color="C1")
xlim([Input1040[1] * Input1550[1] * 0.5, xlim()[2]])
ylim([Output[1] * 0.5, ylim()[2]])
gca()[:set_xscale]("log", nonposx="clip")
gca()[:set_yscale]("log", nonposy="clip")
xlabel("Input power product (W\$^2\$)")
ylabel("Output Power (W)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
