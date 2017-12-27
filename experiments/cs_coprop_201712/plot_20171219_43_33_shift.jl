#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20171219_172640.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20171219_192826.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function gen_selector(site)
    function selector(logicals)
        @assert size(logicals, 2) >= site
        if logicals[1, site] == 0
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[2, site]]
    end
end

data_48 = NaCsData.select_count(params_a, logicals_a, gen_selector(1))
data_56 = NaCsData.select_count(params_b, logicals_b, gen_selector(2))

const spec = OrderedDict(
    :trap155=>linspace(-100, 150, 26),
    :trap50=>linspace(-100, 150, 26),
)

const split_48 = NaCsData.split_data(data_48, spec)
const split_56 = NaCsData.split_data(data_56, spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20171219")

function rabiLine(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end

function gen_model(t)
    return (x, p) -> p[1] .* (1 .- p[2] .* rabiLine.(2π .* (x .- p[3]), t, p[4]))
end

plotx = linspace(-100, 150, 10000)

function fit(t, data, p0, plotx)
    params, ratios, uncs = NaCsData.get_values(data)
    model = gen_model(t)
    fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2], p0)
    return fit.param, estimate_errors(fit), model.(plotx, (fit.param,))
end

fit_48_155 = fit(14.5e-3, split_48[:trap155], [0.9, 0.9, 30, 20], plotx)
# @show Unc.(fit_48_155[1], fit_48_155[2])
fit_48_50 = fit(14.5e-3, split_48[:trap50], [0.9, 0.9, 30, 20], plotx)
# @show Unc.(fit_48_50[1], fit_48_50[2])

fit_56_155 = fit(14.5e-3, split_56[:trap155], [0.9, 0.9, 30, 20], plotx)
# @show Unc.(fit_56_155[1], fit_56_155[2])
fit_56_50 = fit(14.5e-3, split_56[:trap50], [0.9, 0.9, 30, 20], plotx)
# @show Unc.(fit_56_50[1], fit_56_50[2])

figure()
NaCsPlot.plot_survival_data(split_48[:trap155], fmt="C1o", label="15.5mW")
NaCsPlot.plot(plotx, fit_48_155[3], "C1-")
NaCsPlot.plot_survival_data(split_48[:trap50], fmt="C0o", label="5mW")
NaCsPlot.plot(plotx, fit_48_50[3], "C0-")
grid()
xlim([-120, 200])
ylim([0, 1])
text(-115, 0.05, "$(Unc(fit_48_155[1][3], fit_48_155[2][3])) kHz", color="C1")
text(70, 0.4, "$(Unc(fit_48_50[1][3], fit_48_50[2][3])) kHz", color="C0")
legend()
title("Trap at 48MHz")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_48")

figure()
NaCsPlot.plot_survival_data(split_56[:trap155], fmt="C1o", label="15.5mW")
NaCsPlot.plot(plotx, fit_56_155[3], "C1-")
NaCsPlot.plot_survival_data(split_56[:trap50], fmt="C0o", label="5mW")
NaCsPlot.plot(plotx, fit_56_50[3], "C0-")
grid()
xlim([-120, 200])
ylim([0, 1])
text(-115, 0.05, "$(Unc(fit_56_155[1][3], fit_56_155[2][3])) kHz", color="C1")
text(70, 0.4, "$(Unc(fit_56_50[1][3], fit_56_50[2][3])) kHz", color="C0")
legend()
title("Trap at 56.5MHz")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_56")

NaCsPlot.maybe_show()
