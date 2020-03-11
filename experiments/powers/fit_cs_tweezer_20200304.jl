#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci

param_str(fit, i) = @sprintf("%.4f", fit.param[i])
params_strs(fit, suffix="") =
    "\$a$suffix=$(param_str(fit, 1))\$\n\$b$suffix=$(param_str(fit, 2))\$\n\$c$suffix=$(param_str(fit, 3))\$"
params_strs2(fit, suffix="") =
    "\$a$suffix=$(param_str(fit, 1))\$\n\$b$suffix=$(param_str(fit, 2))\$\n\$c$suffix=$(param_str(fit, 3))\$\n\$d$suffix=$(param_str(fit, 4))\$"

const iname_paaom = joinpath(@__DIR__, "data", "cs_tweezer_paaom_20200304.csv")
const data_paaom = readdlm(iname_paaom, ',', Float64, skipstart=1)
const iname_padpaom = joinpath(@__DIR__, "data", "cs_tweezer_padpaom_20200304_345.csv")
const data_padpaom = readdlm(iname_padpaom, ',', Float64, skipstart=1)
const iname_padpaom2 = joinpath(@__DIR__, "data", "cs_tweezer_padpaom2_20200304_345.csv")
const data_padpaom2 = readdlm(iname_padpaom2, ',', Float64, skipstart=1)

@assert data_padpaom2[:, 1] == data_padpaom[:, 1]
data_padpaom2[:, 2] .-= data_padpaom[:, 2]

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

const prefix = joinpath(@__DIR__, "imgs", "cs_tweezer_20200304")

model(x, p) = p[1] .* sin.(p[2] .* sin.(p[3] .* x).^2)
model2(x, p) = p[1] .* sin.(p[2] .* atan.(p[3] .* x).^2).^2
model3(x, p) = p[1] .* (1 .- p[4] .* sin.(p[2] .* atan.(p[3] .* x).^2)).^2

function find_crossing_23(fit2, fit3, maxamp)
    diff23(x) = model2(x, fit2.param) - model3(x, fit3.param)
    a0 = 0.0
    a1 = Float64(maxamp)
    v0 = diff23(a0)
    v1 = diff23(a1)
    assert(v0 * v1 < 0)
    while abs(a0 - a1) > 0.000001
        a2 = (a1 * v0 - a0 * v1) / (v0 - v1)
        v2 = diff23(a2)
        if v2 == 0
            return a2
        end
        a0, a1 = a1, a2
        v0, v1 = v1, v2
    end
    return a1, (model2(a1, fit2.param) + model3(a1, fit3.param)) / 2
end

function find_min_23(fit2, fit3, maxamp)
    sum23(x) = model2(x, fit2.param) + model3(x, fit3.param)
    w1 = (sqrt(5) - 1) / 2
    w2 = 1 - w1

    a0 = 0.0
    a3 = Float64(maxamp)
    a1 = a0 * w1 + a3 * w2
    a2 = a0 * w2 + a3 * w1
    v0 = sum23(a0)
    v1 = sum23(a1)
    v2 = sum23(a2)
    v3 = sum23(a3)
    assert(min(v1, v2) < max(v0, v3))
    while abs(a0 - a1) > 0.00001
        if v1 > v2
            a0, a1, a2, a3 = a1, a2, a1 * w2 + a3 * w1, a3
            v0, v1, v2, v3 = v1, v2, sum23(a2), v3
        else
            a0, a1, a2, a3 = a0, a0 * w1 + a2 * w2, a1, a2
            v0, v1, v2, v3 = v0, sum23(a1), v1, v2
        end
    end
    return a1 > a2 ? (a2, v2) : (a1, v2)
end

plotmax_paaom = maximum(data_paaom[:, 1]) * 1.1
fit_paaom = fit_data(model, data_paaom[:, 1], data_paaom[:, 2],
                     [maximum(data_paaom[:, 2]), 1.5, 2.1], plot_hi=plotmax_paaom)
figure()
plot(data_paaom[:, 1], data_paaom[:, 2], "C2o")
plot(fit_paaom.plotx, fit_paaom.ploty, "C0")
title("PAAOM")
grid()
xlim([0, plotmax_paaom])
ylim([0, ylim()[2]])
text(0.3, 24, "\$a_1\\cdot\\sin(b_1\\cdot\\sin^2(c_1\\cdot AMP))\$", fontsize=18, color="C0")
text(0.45, 2, params_strs(fit_paaom, "_1"), fontsize=20, color="C0")
xlabel("PAAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_paaom")

plotmax_padpaom = maximum(data_padpaom[:, 1]) * 1.05
plotamp_padpaom = linspace(0, plotmax_padpaom, 10000)
fit_padpaom = fit_data(model2, data_padpaom[1:9, 1], data_padpaom[1:9, 2],
                       [maximum(data_padpaom[1:9, 2]), 1.6, 2.1],
                       plotx=plotamp_padpaom)
fit_padpaom2 = fit_data(model3, data_padpaom2[1:9, 1], data_padpaom2[1:9, 2],
                        [maximum(data_padpaom2[1:9, 2]), 1.6, 2.1, 1.0],
                        plotx=plotamp_padpaom)
x0, y0 = find_crossing_23(fit_padpaom, fit_padpaom2, 0.6)
xmin, ymin = find_min_23(fit_padpaom, fit_padpaom2, 0.6)
figure()
plot(data_padpaom[:, 1], data_padpaom[:, 2], "C2o")
plot(fit_padpaom.plotx, fit_padpaom.ploty, "C0-", label="1st order")
plot(data_padpaom2[:, 1], data_padpaom2[:, 2], "C3o")
plot(fit_padpaom2.plotx, fit_padpaom2.ploty, "C1-", label="0th order")
plot(fit_padpaom.plotx, fit_padpaom.ploty .+ fit_padpaom2.ploty,
     "--", color="#d452d4", label="Both")
plot(x0, y0, "bs")
plot(xmin, ymin, "rs")
ax = gca()
ax[:annotate]("x=$(@sprintf("%.5f", x0))\ny=$(@sprintf("%.4f", y0))",
              (x0, y0 + 0.3), xytext=(0.17, 63),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"b"),
              color="b", fontsize="small")
ax[:annotate]("\$x_{min}=$(@sprintf("%.5f", xmin))\$\n\$y_{min}=$(@sprintf("%.4f", ymin))\$",
              (xmin, ymin + 0.3), xytext=(0.55, 63),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"r"),
              color="r", fontsize="small")
title("PADPAOM")
legend(ncol=2, fontsize="small")
grid()
xlim([0, plotmax_padpaom])
ylim([0, 1.3 * max(maximum(data_padpaom[:, 2]), maximum(data_padpaom2[:, 2]))])
text(0.002, 25, "\$a_2\\!\\cdot\\!\\sin^2(b_2\\!\\cdot\\! \\mathrm{atan}^2(c_2\\!\\cdot\\! AMP))\$",
     fontsize=13, color="C0")
text(0.006, 3, params_strs(fit_padpaom, "_2"), fontsize=15, color="C0")
text(0.40, 35, "\$a_3\\!\\cdot\\!(1\\!-\\!d_3\\!\\cdot\\!\\sin(b_3\\!\\cdot\\! \\mathrm{atan}^2(c_3\\!\\cdot\\! AMP)))^2\$", fontsize=13, color="C1")
text(0.66, 6, params_strs2(fit_padpaom2, "_3"), fontsize=15, color="C1")
xlabel("PADPAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_padpaom")

NaCsPlot.maybe_show()
