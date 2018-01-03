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

function fit(model, data, p0, plotx; use_unc=false)
    params, ratios, uncs = NaCsData.get_values(data)
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2], p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
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
    NaCsPlot.plot_survival_data(split_50_mf3_lo[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_50_mf3_hi[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_50_mf4_lo[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_50_mf4_hi[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_62_mf3_lo[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_62_mf3_hi[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_62_mf4_lo[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
    NaCsPlot.plot_survival_data(split_62_mf4_hi[i], yoffset=yoffset, fmt="C$(i - 1)o")
    plot(plotx, fit_res[3] .+ yoffset, "C$(i - 1)-")
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
bshifts_plot = linspace(-0.35, 0.35, 10000)

const shift_50_mf3 = (res_50_mf3_hi .- res_50_mf3_lo) ./ (15.5 - 5) .* 15.5
const shift_50_mf3_unc = sqrt.(res_50_mf3_hi_unc.^2 .+ res_50_mf3_lo_unc.^2) ./ (15.5 - 5) .* 15.5
const shift_50_mf4 = (res_50_mf4_hi .- res_50_mf4_lo) ./ (15.5 - 5) .* 15.5
const shift_50_mf4_unc = sqrt.(res_50_mf4_hi_unc.^2 .+ res_50_mf4_lo_unc.^2) ./ (15.5 - 5) .* 15.5
const shift_62_mf3 = (res_62_mf3_hi .- res_62_mf3_lo) ./ (15.5 - 5) .* 15.5
const shift_62_mf3_unc = sqrt.(res_62_mf3_hi_unc.^2 .+ res_62_mf3_lo_unc.^2) ./ (15.5 - 5) .* 15.5
const shift_62_mf4 = (res_62_mf4_hi .- res_62_mf4_lo) ./ (15.5 - 5) .* 15.5
const shift_62_mf4_unc = sqrt.(res_62_mf4_hi_unc.^2 .+ res_62_mf4_lo_unc.^2) ./ (15.5 - 5) .* 15.5

function _shift_model(b, idx, p)
    slope = p[5]
    offset = p[idx]
    if idx == 1 || idx == 3
        slope = slope * 6 / 7
    end
    return (b - offset) * slope
end
function shift_model(x::Integer, p)
    nbs = length(BShifts)
    if x > 10
        x += 1
    elseif x >= 23
        x += 2
    end
    return _shift_model(BShifts[(x - 1) % nbs + 1], (x - 1) ÷ nbs + 1, p)
end
shift_model(x, p) = shift_model.(x, (p,))

shift_fit = curve_fit(shift_model, 1:(length(BShifts) * 4 - 2),
                      [shift_50_mf3; shift_50_mf4[[1:3; 5:7]];
                       shift_62_mf3; shift_62_mf4[[1:3; 5:7]]],
                      [0.6, 0.6, 0.6, 0.6, 35])

shift_param = shift_fit.param
shift_param_unc = Unc.(shift_param, estimate_errors(shift_fit))

# shift_50_mf3_param = shift_50_mf3_fit.param
# shift_50_mf3_param_unc = Unc.(shift_50_mf3_param, estimate_errors(shift_50_mf3_fit))
# shift_50_mf4_param = shift_50_mf4_fit.param
# shift_50_mf4_param_unc = Unc.(shift_50_mf4_param, estimate_errors(shift_50_mf4_fit))
# shift_62_mf3_param = shift_62_mf3_fit.param
# shift_62_mf3_param_unc = Unc.(shift_62_mf3_param, estimate_errors(shift_62_mf3_fit))
# shift_62_mf4_param = shift_62_mf4_fit.param
# shift_62_mf4_param_unc = Unc.(shift_62_mf4_param, estimate_errors(shift_62_mf4_fit))

figure()
p = errorbar(BShifts, shift_50_mf3, shift_50_mf3_unc, fmt="o", label="\$50MHz\\ m_F=3\$")
c = p[1][:get_color]()
plot(bshifts_plot, _shift_model.(bshifts_plot, 1, (shift_param,)), "-", color=c)
gcf()[:text](0.75, 0.45,
             "\$(B_y - $(shift_param_unc[1])) * $(shift_param_unc[5] * (6 / 7))\$",
             color=c, fontsize=20, horizontalalignment="left", verticalalignment="center",
             clip_on=false)
p = errorbar(BShifts, shift_50_mf4, shift_50_mf4_unc, fmt="o", label="\$50MHz\\ m_F=4\$")
c = p[1][:get_color]()
plot(bshifts_plot, _shift_model.(bshifts_plot, 2, (shift_param,)), "-", color=c)
gcf()[:text](0.75, 0.35,
             "\$(B_y - $(shift_param_unc[2])) * $(shift_param_unc[5])\$",
             color=c, fontsize=20, horizontalalignment="left", verticalalignment="center",
             clip_on=false)
p = errorbar(BShifts, shift_62_mf3, shift_62_mf3_unc, fmt="o", label="\$62MHz\\ m_F=3\$")
c = p[1][:get_color]()
plot(bshifts_plot, _shift_model.(bshifts_plot, 3, (shift_param,)), "-", color=c)
gcf()[:text](0.75, 0.25,
             "\$(B_y - $(shift_param_unc[3])) * $(shift_param_unc[5] * (6 / 7))\$",
             color=c, fontsize=20, horizontalalignment="left", verticalalignment="center",
             clip_on=false)
p = errorbar(BShifts, shift_62_mf4, shift_62_mf4_unc, fmt="o", label="\$62MHz\\ m_F=4\$")
c = p[1][:get_color]()
plot(bshifts_plot, _shift_model.(bshifts_plot, 4, (shift_param,)), "-", color=c)
gcf()[:text](0.75, 0.15,
             "\$(B_y - $(shift_param_unc[4])) * $(shift_param_unc[5])\$",
             color=c, fontsize=20, horizontalalignment="left", verticalalignment="center",
             clip_on=false)
title("Shift from trap light")
xlabel("\$B_y (G)\$")
ylabel("Shift (kHz)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
grid()
NaCsPlot.maybe_save("$(prefix)_trap_shifts")

const shift_50_diff = shift_50_mf4 .- shift_50_mf3 .* (7 / 6)
const shift_50_diff_unc = sqrt.(shift_50_mf4_unc.^2 .+ shift_50_mf3_unc.^2 .* (7 / 6)^2)
const shift_62_diff = shift_62_mf4 .- shift_62_mf3 .* (7 / 6)
const shift_62_diff_unc = sqrt.(shift_62_mf4_unc.^2 .+ shift_62_mf3_unc.^2 .* (7 / 6)^2)

figure()
p = errorbar(BShifts, shift_50_diff, shift_50_diff_unc, fmt="o-", label="\$50MHz\$")
p = errorbar(BShifts, shift_62_diff, shift_62_diff_unc, fmt="o-", label="\$62MHz\$")
title("Extra shift in \$m_F=4\$")
xlabel("\$B_y (G)\$")
ylabel("Shift (kHz)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
grid()
NaCsPlot.maybe_save("$(prefix)_trap_shifts_extra")

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
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
grid()
NaCsPlot.maybe_save("$(prefix)_zeeman_shifts")

const zeeman_mf3 =
    (zeeman_50_mf3 ./ zeeman_50_mf3_unc.^2 .+ zeeman_62_mf3 ./ zeeman_62_mf3_unc.^2) ./
    (1 ./ zeeman_50_mf3_unc.^2 .+ 1 ./ zeeman_62_mf3_unc.^2)
const zeeman_mf3_unc = 1 ./ sqrt.(1 ./ zeeman_50_mf3_unc.^2 .+ 1 ./ zeeman_62_mf3_unc.^2)
const zeeman_mf4 =
    (zeeman_50_mf4 ./ zeeman_50_mf4_unc.^2 .+ zeeman_62_mf4 ./ zeeman_62_mf4_unc.^2) ./
    (1 ./ zeeman_50_mf4_unc.^2 .+ 1 ./ zeeman_62_mf4_unc.^2)
const zeeman_mf4_unc = 1 ./ sqrt.(1 ./ zeeman_50_mf4_unc.^2 .+ 1 ./ zeeman_62_mf4_unc.^2)

function _zeemans_model(b, ismf3, p)
    if ismf3
        δgF = 350 * 6
        offset = p[3]
    else
        δgF = 350 * 7
        offset = p[4]
    end
    B0 = 8.8392
    return sqrt(((p[1] * b - p[2]))^2 + B0^2) * δgF - B0 * δgF + offset
end
function zeemans_model(x::Integer, p)
    if x <= length(BShifts)
        return _zeemans_model(BShifts[x], true, p)
    else
        return _zeemans_model(BShifts[x - 7], false, p)
    end
end
zeemans_model(x, p) = zeemans_model.(x, (p,))

zeemans_fit = curve_fit(zeemans_model, 1:(length(BShifts) * 2),
                        [zeeman_mf3; zeeman_mf4], [1 ./ zeeman_mf3_unc; 1 ./ zeeman_mf4_unc],
                        [1.0, 0.0, 0.0, 0.0])
zeemans_param = zeemans_fit.param
zeemans_unc = estimate_errors(zeemans_fit)
@show Unc.(zeemans_param, zeemans_unc)

figure()
errorbar(BShifts, zeeman_mf3, zeeman_mf3_unc, fmt="C0o", label="\$m_F=3\$")
plot(bshifts_plot, _zeemans_model.(bshifts_plot, true, (zeemans_fit.param,)), "C0-")
errorbar(BShifts, zeeman_mf4, zeeman_mf4_unc, fmt="C1o", label="\$m_F=4\$")
plot(bshifts_plot, _zeemans_model.(bshifts_plot, false, (zeemans_fit.param,)), "C1-")
gcf()[:text](0.83, 0.35,
             "\$B_0=$(Unc(zeemans_param[2] * 1000, zeemans_unc[2] * 1000, Sci))\\ mG\$",
             color="k", fontsize=20, horizontalalignment="left", verticalalignment="center",
             clip_on=false)
gcf()[:text](0.78, 0.2,
             "\$\\dfrac{B_{real}}{B_z}=$(Unc(zeemans_param[1], zeemans_unc[1], Sci))\$",
             color="k", fontsize=20, horizontalalignment="left", verticalalignment="center",
             clip_on=false)
title("Shift from B field")
xlabel("\$B_y (G)\$")
ylabel("Shift (kHz)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
grid()
NaCsPlot.maybe_save("$(prefix)_zeeman_shifts_avg")

NaCsPlot.maybe_show()
