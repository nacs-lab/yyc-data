#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180121_145134.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :x=>(linspace(18.11, 18.20, 10), linspace(19.05, 19.20, 16),
         linspace(19.54, 19.69, 16)),
    :y=>(linspace(18.10, 18.20, 11), linspace(19.04, 19.20, 17),
         linspace(19.52, 19.71, 20)),
    :z=>(linspace(18.49, 18.545, 12), linspace(18.57, 18.635, 14),
         linspace(18.66, 18.71, 11), linspace(18.75, 18.80, 11),
         linspace(18.83, 18.89, 13), linspace(18.91, 18.97, 13),
         linspace(18.99, 19.06, 15)),
    :x0=>linspace(18.50, 18.70, 11),
    :y0=>linspace(18.50, 18.70, 11),
)

function fit(model, data, p0; plotx=nothing, use_unc=false, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        plotx = linspace(mid - span * plot_scale / 2, mid + span * plot_scale / 2, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2].^2, p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

model(x, p) = p[1] .* exp.(.-((x .- p[2]) ./ p[3]).^2)

const split_a = NaCsData.split_data(data_a, spec)

const data_x = split_a[:x]
const data_y = split_a[:y]
const data_z = split_a[:z]
const data_x0 = split_a[:x0]
const data_y0 = split_a[:y0]

function compute_pg(fit_p1, fit_m1)
    hp1 = fit_p1.param[1]
    up1 = fit_p1.unc[1]
    hm1 = fit_m1.param[1]
    um1 = fit_m1.unc[1]
    r = hm1 / hp1
    ur = r * sqrt((up1 / hp1)^2 + (um1 / hm1)^2)
    pg = 1 - r
    upg = ur
    return Unc(pg, upg)
end

const prefix = joinpath(@__DIR__, "imgs", "data_20180121_145134_spectrum_cold")

fit_xp1 = fit(model, data_x[1][4:end - 2], [0.88, 18.162, 0.02], plot_scale=2)
@show Unc.(fit_xp1.param, fit_xp1.unc)
fit_xm1 = fit(model, data_x[2][6:end - 2], [0.02, 19.140, 0.02], plot_scale=2)
@show Unc.(fit_xm1.param, fit_xm1.unc)
pg_x = compute_pg(fit_xp1, fit_xm1)
figure()
NaCsPlot.plot_survival_data(data_x[1], fmt="C0.")
plot(fit_xp1.plotx, fit_xp1.ploty, color="C0")
NaCsPlot.plot_survival_data(data_x[2], fmt="C0.")
plot(fit_xm1.plotx, fit_xm1.ploty, color="C0")
NaCsPlot.plot_survival_data(data_x[3], fmt="C0.-")
text(18.50, 0.6, "\$P_{ground}=$(pg_x * 100) \\%\$")
grid()
ylim([0, 1])
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

fit_yp1 = fit(model, data_y[1][4:end - 2], [0.89, 18.140, 0.02], plot_scale=2.5)
@show Unc.(fit_yp1.param, fit_yp1.unc)
fit_ym1 = fit(model, data_y[2][6:end - 2], [0.02, 19.140, 0.03], plot_scale=2)
@show Unc.(fit_ym1.param, fit_ym1.unc)
pg_y = compute_pg(fit_yp1, fit_ym1)
figure()
NaCsPlot.plot_survival_data(data_y[1], fmt="C0.")
plot(fit_yp1.plotx, fit_yp1.ploty, color="C0")
NaCsPlot.plot_survival_data(data_y[2], fmt="C0.")
plot(fit_ym1.plotx, fit_ym1.ploty, color="C0")
NaCsPlot.plot_survival_data(data_y[3], fmt="C0.-")
text(18.50, 0.6, "\$P_{ground}=$(pg_y * 100) \\%\$")
grid()
ylim([0, 1])
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

fit_zp1 = fit(model, data_z[1][3:end - 2], [0.88, 18.515, 0.01], plot_scale=2)
@show Unc.(fit_zp1.param, fit_zp1.unc)
fit_zm1 = fit(model, data_z[3], [0.02, 18.690, 0.01], plot_scale=2)
@show Unc.(fit_zm1.param, fit_zm1.unc)
pg_z = compute_pg(fit_zp1, fit_zm1)
figure()
NaCsPlot.plot_survival_data(data_z[1], fmt="C0.")
plot(fit_zp1.plotx, fit_zp1.ploty, color="C0")
NaCsPlot.plot_survival_data(data_z[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[3], fmt="C0.")
plot(fit_zm1.plotx, fit_zm1.ploty, color="C0")
NaCsPlot.plot_survival_data(data_z[4], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[5], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[6], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[7], fmt="C0.-")
text(18.60, 0.7, "\$P_{ground}=$(pg_z * 100) \\%\$")
grid()
ylim([0, 1])
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

pg_3d = pg_x * pg_y * pg_z
@show pg_3d

figure()
NaCsPlot.plot_survival_data(data_x0, fmt="C0.-")
grid()
ylim([0, 1])
title("X carrier with square pulse")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x0")

figure()
NaCsPlot.plot_survival_data(data_y0, fmt="C0.-")
grid()
ylim([0, 1])
title("Y carrier with square pulse")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y0")

NaCsPlot.maybe_show()
