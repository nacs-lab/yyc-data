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
    "\$a$suffix=$(param_str(fit, 1))\$\n\$b$suffix=$(param_str(fit, 2))\$"

const iname_paaom = joinpath(@__DIR__, "data", "cs_tweezer_paaom_20200207.csv")
const data_paaom = readdlm(iname_paaom, ',', Float64, skipstart=1)
const iname_padpaom = joinpath(@__DIR__, "data", "cs_tweezer_padpaom_20200207.csv")
const data_padpaom = readdlm(iname_padpaom, ',', Float64, skipstart=1)
const iname_padpaom2 = joinpath(@__DIR__, "data", "cs_tweezer_padpaom2_20200207.csv")
const data_padpaom2 = readdlm(iname_padpaom2, ',', Float64, skipstart=1)

@assert data_padpaom2[:, 1] == data_padpaom[:, 1]
data_padpaom2[:, 2] .-= data_padpaom[:, 2]

const prefix = joinpath(@__DIR__, "imgs", "cs_tweezer_20200207")

model(x, p) = p[1] .* sin.(p[2] .* sin.(p[3] .* x).^2)
model2(x, p) = p[1] .* sin.(p[2] .* x).^4
model3(x, p) = p[1] .* cos.(p[2] .* x).^4
model2′(x, p) = (4 * p[1] * p[2]) .* sin.(p[2] .* x).^3 .* cos.(p[2] .* x)
model3′(x, p) = (-4 * p[1] * p[2]) .* cos.(p[2] .* x).^3 .* sin.(p[2] .* x)

model2′′(x, p) = (4 * p[1] * p[2]^2) .*
    (3 .* sin.(p[2] .* x).^2 .* cos.(p[2] .* x).^2 .- sin.(p[2] .* x).^4)
model3′′(x, p) = (-4 * p[1] * p[2]^2) .* (-3 .* cos.(p[2] .* x).^2 .* sin.(p[2] .* x).^2 .+
                                            cos.(p[2] .* x).^4)

model23′(x, p2, p3) = model2′(x, p2) .+ model3′(x, p3)
model23′′(x, p2, p3) = model2′′(x, p2) .+ model3′′(x, p3)

diff23(x, p2, p3) = model2(x, p2) .- model3(x, p3)
diff23′(x, p2, p3) = model2′(x, p2) .- model3′(x, p3)

function solve23(x0, p2, p3)
    while true
        y = diff23(x0, p2, p3)
        y′ = diff23′(x0, p2, p3)
        dx = y / y′
        x1 = x0 - dx
        if abs(dx) < 1e-8
            return (x1, (model2(x1, p2) + model3(x1, p3)) / 2)
        end
        x0 = x1
    end
end

function solve23′(x0, p2, p3)
    while true
        y = model23′(x0, p2, p3)
        y′ = model23′′(x0, p2, p3)
        dx = y / y′
        x1 = x0 - dx
        if abs(dx) < 1e-8
            return (x1, model2(x1, p2) + model3(x1, p3))
        end
        x0 = x1
    end
end

plotmax_paaom = maximum(data_paaom[:, 1]) * 1.1
plotamp_paaom = linspace(0, plotmax_paaom, 10000)
fit_paaom = curve_fit(model, data_paaom[:, 1], data_paaom[:, 2],
                      [maximum(data_paaom[:, 2]), 1.5, 2.1])
figure()
plot(data_paaom[:, 1], data_paaom[:, 2], "C2o")
plot(plotamp_paaom, model(plotamp_paaom, fit_paaom.param), "C0")
title("PAAOM")
grid()
xlim([0, plotmax_paaom])
ylim([0, ylim()[2]])
text(0.3, 18, "\$a_1\\cdot\\sin(b_1\\cdot\\sin^2(c_1\\cdot AMP))\$", fontsize=18, color="C0")
text(0.45, 1, params_strs(fit_paaom, "_1"), fontsize=20, color="C0")
xlabel("PAAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_paaom")

plotmax_padpaom = maximum(data_padpaom[:, 1]) * 1.05
plotamp_padpaom = linspace(0, plotmax_padpaom, 10000)

fit_padpaom = curve_fit(model2, data_padpaom[:, 1], data_padpaom[:, 2],
                        [maximum(data_padpaom[:, 2]), 1.6])
fit_padpaom2 = curve_fit(model3, data_padpaom2[:, 1], data_padpaom2[:, 2],
                         [maximum(data_padpaom2[:, 2]), 1.0])
x0, y0 = solve23(0.3, fit_padpaom.param, fit_padpaom2.param)
xmin, ymin = solve23′(0.3, fit_padpaom.param, fit_padpaom2.param)
figure()
plot(data_padpaom[:, 1], data_padpaom[:, 2], "C2o")
plot(plotamp_padpaom, model2(plotamp_padpaom, fit_padpaom.param), "C0-", label="1st order")
plot(data_padpaom2[:, 1], data_padpaom2[:, 2], "C3o")
plot(plotamp_padpaom, model3(plotamp_padpaom, fit_padpaom2.param), "C1-", label="0th order")
plot(plotamp_padpaom, model2(plotamp_padpaom, fit_padpaom.param) .+
     model3(plotamp_padpaom, fit_padpaom2.param), "--", color="#d452d4", label="Both")
plot(x0, y0, "bs")
plot(xmin, ymin, "rs")
ax = gca()
ax[:annotate]("x=$(@sprintf("%.5f", x0))\ny=$(@sprintf("%.4f", y0))",
              (x0, y0 + 0.3), xytext=(0.64, 28),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"b"),
              color="b")
ax[:annotate]("\$x_{min}=$(@sprintf("%.5f", xmin))\$\n\$y_{min}=$(@sprintf("%.4f", ymin))\$",
              (xmin, ymin + 0.3), xytext=(0.21, 45),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"r"),
              color="r")
title("PADPAOM")
legend(ncol=2, fontsize="small")
grid()
xlim([0, plotmax_padpaom])
ylim([0, 1.3 * max(maximum(data_padpaom[:, 2]), maximum(data_padpaom2[:, 2]))])
text(0.002, 18, "\$a_2\\cdot\\sin^4(b_2\\cdot AMP)\$", fontsize=18, color="C0")
text(0.006, 3, params_strs2(fit_padpaom, "_2"), fontsize=20, color="C0")
text(0.585, 18, "\$a_3\\cdot\\cos^4(b_3\\cdot AMP)\$", fontsize=18, color="C1")
text(0.62, 3, params_strs2(fit_padpaom2, "_3"), fontsize=20, color="C1")
xlabel("PADPAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_padpaom")

NaCsPlot.maybe_show()
