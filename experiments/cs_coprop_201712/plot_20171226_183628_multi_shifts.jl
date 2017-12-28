#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20171226_183628.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function gen_selector(site)
    function selector(logicals)
        @assert size(logicals, 2) >= site
        if logicals[1, site] == 0
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[2, site]]
    end
end

data_50 = NaCsData.select_count(params_a, logicals_a, gen_selector(1))
data_62 = NaCsData.select_count(params_a, logicals_a, gen_selector(2))

data_50_mf3 = data_50[1:2:end]
data_50_mf4 = data_50[2:2:end]

data_62_mf3 = data_62[1:2:end]
data_62_mf4 = data_62[2:2:end]

data_50_mf3_lo = data_50_mf3[1:end ÷ 2]
data_50_mf3_hi = data_50_mf3[end ÷ 2 + 1:end]
data_50_mf4_lo = data_50_mf4[1:end ÷ 2]
data_50_mf4_hi = data_50_mf4[end ÷ 2 + 1:end]

data_62_mf3_lo = data_62_mf3[1:end ÷ 2]
data_62_mf3_hi = data_62_mf3[end ÷ 2 + 1:end]
data_62_mf4_lo = data_62_mf4[1:end ÷ 2]
data_62_mf4_hi = data_62_mf4[end ÷ 2 + 1:end]

const spec = ((linspace(-80, 100, 19) for _ in 1:7)...)

const split_50_mf3_lo = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_50_mf3_lo), spec)
const split_50_mf3_hi = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_50_mf3_hi), spec)
const split_50_mf4_lo = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_50_mf4_lo), spec)
const split_50_mf4_hi = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_50_mf4_hi), spec)
const split_62_mf3_lo = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_62_mf3_lo), spec)
const split_62_mf3_hi = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_62_mf3_hi), spec)
const split_62_mf4_lo = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_62_mf4_lo), spec)
const split_62_mf4_hi = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_62_mf4_hi), spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20171226")

function rabiLine(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end

function gen_model3(t)
    return (x, p) -> p[1] .* (1 .- p[2] .* rabiLine.(2π .* (x .- p[3]), t, p[4]))
end
function gen_model4(t)
    return (x, p) -> p[1] .* rabiLine.(2π .* (x .- p[2]), t, p[3])
end

const plotx = linspace(-90, 110, 10000)

function fit(model, data, p0, plotx)
    params, ratios, uncs = NaCsData.get_values(data)
    # fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2], p0)
    fit = curve_fit(model, params, ratios[:, 2], p0)
    return fit.param, estimate_errors(fit), model.(plotx, (fit.param,))
end

const res_50_mf3_lo = Float64[]
const res_50_mf3_lo_unc = Float64[]
const res_50_mf3_hi = Float64[]
const res_50_mf3_hi_unc = Float64[]
const res_50_mf4_lo = Float64[]
const res_50_mf4_lo_unc = Float64[]
const res_50_mf4_hi = Float64[]
const res_50_mf4_hi_unc = Float64[]
const res_62_mf3_lo = Float64[]
const res_62_mf3_lo_unc = Float64[]
const res_62_mf3_hi = Float64[]
const res_62_mf3_hi_unc = Float64[]
const res_62_mf4_lo = Float64[]
const res_62_mf4_lo_unc = Float64[]
const res_62_mf4_hi = Float64[]
const res_62_mf4_hi_unc = Float64[]

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1) - 0.15
    fit_res = fit(gen_model3(15e-3), split_50_mf3_lo[i], [0.8, 0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_50_mf3_lo[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_50_mf3_lo, fit_res[1][end - 1])
    push!(res_50_mf3_lo_unc, fit_res[2][end - 1])
end
title("\$50MHz\\  5mW\\  m_F=3\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2])
NaCsPlot.maybe_save("$(prefix)_50_mf3_lo")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1) - 0.15
    fit_res = fit(gen_model3(15e-3), split_50_mf3_hi[i], [0.8, 0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_50_mf3_hi[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_50_mf3_hi, fit_res[1][end - 1])
    push!(res_50_mf3_hi_unc, fit_res[2][end - 1])
end
title("\$50MHz\\  15.5mW\\  m_F=3\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2])
NaCsPlot.maybe_save("$(prefix)_50_mf3_hi")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1)
    fit_res = fit(gen_model4(15e-3), split_50_mf4_lo[i], [0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_50_mf4_lo[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_50_mf4_lo, fit_res[1][end - 1])
    push!(res_50_mf4_lo_unc, fit_res[2][end - 1])
end
title("\$50MHz\\ 5mW\\  m_F=4\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2.2])
NaCsPlot.maybe_save("$(prefix)_50_mf4_lo")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1)
    fit_res = fit(gen_model4(15e-3), split_50_mf4_hi[i], [0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_50_mf4_hi[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_50_mf4_hi, fit_res[1][end - 1])
    push!(res_50_mf4_hi_unc, fit_res[2][end - 1])
end
title("\$50MHz\\ 15.5mW\\  m_F=4\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2.2])
NaCsPlot.maybe_save("$(prefix)_50_mf4_hi")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1) - 0.15
    fit_res = fit(gen_model3(15e-3), split_62_mf3_lo[i], [0.8, 0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_62_mf3_lo[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_62_mf3_lo, fit_res[1][end - 1])
    push!(res_62_mf3_lo_unc, fit_res[2][end - 1])
end
title("\$62MHz\\  5mW\\  m_F=3\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2])
NaCsPlot.maybe_save("$(prefix)_62_mf3_lo")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1) - 0.15
    fit_res = fit(gen_model3(15e-3), split_62_mf3_hi[i], [0.8, 0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_62_mf3_hi[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_62_mf3_hi, fit_res[1][end - 1])
    push!(res_62_mf3_hi_unc, fit_res[2][end - 1])
end
title("\$62MHz\\  15.5mW\\  m_F=3\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2])
NaCsPlot.maybe_save("$(prefix)_62_mf3_hi")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1)
    fit_res = fit(gen_model4(15e-3), split_62_mf4_lo[i], [0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_62_mf4_lo[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_62_mf4_lo, fit_res[1][end - 1])
    push!(res_62_mf4_lo_unc, fit_res[2][end - 1])
end
title("\$62MHz\\ 5mW\\  m_F=4\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2.2])
NaCsPlot.maybe_save("$(prefix)_62_mf4_lo")

figure()
for i in 1:7
    yoffset = 0.2 * (i - 1)
    fit_res = fit(gen_model4(15e-3), split_62_mf4_hi[i], [0.8, 30, 20], plotx)
    NaCsPlot.plot_survival_data(split_62_mf4_hi[i], yoffset=yoffset, fmt="C$(i)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i)-")
    push!(res_62_mf4_hi, fit_res[1][end - 1])
    push!(res_62_mf4_hi_unc, fit_res[2][end - 1])
end
title("\$62MHz\\ 15.5mW\\  m_F=4\$")
xlabel("Detuning (kHz)")
ylabel("Survival (shifted)")
grid()
ylim([0, 2.2])
NaCsPlot.maybe_save("$(prefix)_62_mf4_hi")

const BShifts = linspace(-0.3, 0.3, 7)

const shift_50_mf3 = (res_50_mf3_hi .- res_50_mf3_lo) ./ (15.5 - 5) .* 15.5
const shift_50_mf3_unc = sqrt.(res_50_mf3_hi_unc.^2 .+ res_50_mf3_lo_unc.^2) ./ (15.5 - 5) .* 15.5
const shift_50_mf4 = (res_50_mf4_hi .- res_50_mf4_lo) ./ (15.5 - 5) .* 15.5
const shift_50_mf4_unc = sqrt.(res_50_mf4_hi_unc.^2 .+ res_50_mf4_lo_unc.^2) ./ (15.5 - 5) .* 15.5
const shift_62_mf3 = (res_62_mf3_hi .- res_62_mf3_lo) ./ (15.5 - 5) .* 15.5
const shift_62_mf3_unc = sqrt.(res_62_mf3_hi_unc.^2 .+ res_62_mf3_lo_unc.^2) ./ (15.5 - 5) .* 15.5
const shift_62_mf4 = (res_62_mf4_hi .- res_62_mf4_lo) ./ (15.5 - 5) .* 15.5
const shift_62_mf4_unc = sqrt.(res_62_mf4_hi_unc.^2 .+ res_62_mf4_lo_unc.^2) ./ (15.5 - 5) .* 15.5

figure()
errorbar(BShifts, shift_50_mf3, shift_50_mf3_unc, label="\$50MHz\\ m_F=3\$")
errorbar(BShifts, shift_50_mf4, shift_50_mf4_unc, label="\$50MHz\\ m_F=4\$")
errorbar(BShifts, shift_62_mf3, shift_62_mf3_unc, label="\$62MHz\\ m_F=3\$")
errorbar(BShifts, shift_62_mf4, shift_62_mf4_unc, label="\$62MHz\\ m_F=4\$")
title("Shift from trap light")
xlabel("\$B_y (G)\$")
ylabel("Shift (kHz)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
grid()
NaCsPlot.maybe_save("$(prefix)_trap_shifts")

const zeeman_50_mf3 = res_50_mf3_lo .* (15.5 / (15.5 - 5)) .+ res_50_mf3_hi .* (-5 / (15.5 - 5))
const zeeman_50_mf3_unc = sqrt.((res_50_mf3_lo_unc .* (15.5 / (15.5 - 5))).^2 .+
                                (res_50_mf3_hi_unc .* (-5 / (15.5 - 5))).^2)
const zeeman_50_mf4 = res_50_mf4_lo .* (15.5 / (15.5 - 5)) .+ res_50_mf4_hi .* (-5 / (15.5 - 5))
const zeeman_50_mf4_unc = sqrt.((res_50_mf4_lo_unc .* (15.5 / (15.5 - 5))).^2 .+
                                (res_50_mf4_hi_unc .* (-5 / (15.5 - 5))).^2)
const zeeman_62_mf3 = res_62_mf3_lo .* (15.5 / (15.5 - 5)) .+ res_62_mf3_hi .* (-5 / (15.5 - 5))
const zeeman_62_mf3_unc = sqrt.((res_62_mf3_lo_unc .* (15.5 / (15.5 - 5))).^2 .+
                                (res_62_mf3_hi_unc .* (-5 / (15.5 - 5))).^2)
const zeeman_62_mf4 = res_62_mf4_lo .* (15.5 / (15.5 - 5)) .+ res_62_mf4_hi .* (-5 / (15.5 - 5))
const zeeman_62_mf4_unc = sqrt.((res_62_mf4_lo_unc .* (15.5 / (15.5 - 5))).^2 .+
                                (res_62_mf4_hi_unc .* (-5 / (15.5 - 5))).^2)

figure()
errorbar(BShifts, zeeman_50_mf3, zeeman_50_mf3_unc, label="\$50MHz\\ m_F=3\$")
errorbar(BShifts, zeeman_50_mf4, zeeman_50_mf4_unc, label="\$50MHz\\ m_F=4\$")
errorbar(BShifts, zeeman_62_mf3, zeeman_62_mf3_unc, label="\$62MHz\\ m_F=3\$")
errorbar(BShifts, zeeman_62_mf4, zeeman_62_mf4_unc, label="\$62MHz\\ m_F=4\$")
title("Shift from B field")
xlabel("\$B_y (G)\$")
ylabel("Shift (kHz)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
grid()
NaCsPlot.maybe_save("$(prefix)_zeeman_shifts")

NaCsPlot.maybe_show()
