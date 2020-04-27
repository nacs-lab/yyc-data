#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200406_234752.mat",
                "data_20200407_161453.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [([0.0, 2, 4, 10, 20],
                [0.0, 20, 40, 80, 160],
                [0.0, 20, 40, 80, 160] .* 2,
                [0.0, 20, 50],
                [0.0, 80, 160],
                [0.0, 80, 160] .* 2),
               ([0.0, 4, 10],
                [0.0, 40, 80],
                [0.0, 80, 160] .* 2,
                [0.0, 50],
                [0.0, 160],
                [0.0, 160] .* 2)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4,)), maxcnts, specs)

const data_15a = [datas_nacs[1][4]; datas_nacs[2][4]]
const data_15m = [datas_nacs[1][1]; datas_nacs[2][1]]
const data_6a = [datas_nacs[1][5]; datas_nacs[2][5]]
const data_6m = [datas_nacs[1][2]; datas_nacs[2][2]]
const data_3a = [datas_nacs[1][6]; datas_nacs[2][6]]
const data_3m = [datas_nacs[1][3]; datas_nacs[2][3]]

const data_15cs = [datas_cs[1][4]; datas_cs[2][4]]
const data_6cs = [datas_cs[1][5]; datas_cs[2][5]]
const data_3cs = [datas_cs[1][6]; datas_cs[2][6]]
const data_15na = [datas_na[1][4]; datas_na[2][4]; datas_na[1][1]; datas_na[2][1]]
const data_6na = [datas_na[1][5]; datas_na[2][5]; datas_na[1][2]; datas_na[2][2]]
const data_3na = [datas_na[1][6]; datas_na[2][6]; datas_na[1][3]; datas_na[2][3]]

const prefix = joinpath(@__DIR__, "imgs", "data_20200406_234752_pa_time_3322")

function model_expoff(x, p)
    p[3] .+ p[1] .* exp.(.- x ./ p[2])
end
function model_exp(x, p)
    p[1] .* exp.(.- x .* p[2])
end
function gen_power_model(pwr)
    return function (x, p)
        return p[1] .* x.^pwr
    end
end

const fit_15a = fit_survival(model_exp, data_15a, [0.3, 0.05])
const fit_15m = fit_survival(model_exp, data_15m, [0.3, 0.05])
const fit_6a = fit_survival(model_exp, data_6a, [0.3, 0.05])
const fit_6m = fit_survival(model_exp, data_6m, [0.3, 0.05])
const fit_3a = fit_survival(model_exp, data_3a, [0.3, 0.05])
const fit_3m = fit_survival(model_exp, data_3m, [0.3, 0.05])

const fit_15cs = fit_survival(model_exp, data_15cs, [0.8, 0.05])
const fit_6cs = fit_survival(model_exp, data_6cs, [0.8, 0.05])
const fit_3cs = fit_survival(model_exp, data_3cs, [0.8, 0.05])

const fit_15na = fit_survival(model_exp, data_15na, [0.8, 0.05])
const fit_6na = fit_survival(model_exp, data_6na, [0.8, 0.05])
const fit_3na = fit_survival(model_exp, data_3na, [0.8, 0.05])

# const fit_15a2 = fit_survival(model_exp, data_15a2, [0.3, 0.05])
# const fit_15m2 = fit_survival(model_exp, data_15m2, [0.3, 0.05])
# const fit_6a2 = fit_survival(model_exp, data_6a2, [0.3, 0.05])
# const fit_6m2 = fit_survival(model_exp, data_6m2, [0.3, 0.05])
# const fit_3a2 = fit_survival(model_exp, data_3a2, [0.3, 0.05])
# const fit_3m2 = fit_survival(model_exp, data_3m2, [0.3, 0.05])

# const fit_15cs2 = fit_survival(model_exp, data_15cs2, [0.8, 0.05])
# const fit_6cs2 = fit_survival(model_exp, data_6cs2, [0.8, 0.05])
# const fit_3cs2 = fit_survival(model_exp, data_3cs2, [0.8, 0.05])

# const fit_15na2 = fit_survival(model_exp, data_15na2, [0.8, 0.05])
# const fit_6na2 = fit_survival(model_exp, data_6na2, [0.8, 0.05])
# const fit_3na2 = fit_survival(model_exp, data_3na2, [0.8, 0.05])

figure(figsize=[12.6, 11.2])

const γ_15pa = fit_15m.uncs[2] - fit_15na.uncs[2] - fit_15cs.uncs[2]
subplot(2, 2, 1)
NaCsPlot.plot_survival_data(data_15na, fmt="C0s", label="Na")
plot(fit_15na.plotx, fit_15na.ploty, "C0-")
NaCsPlot.plot_survival_data(data_15cs, fmt="C1s", label="Cs")
plot(fit_15cs.plotx, fit_15cs.ploty, "C1-")
NaCsPlot.plot_survival_data(data_15a, fmt="C2s", label="Hot")
plot(fit_15a.plotx, fit_15a.ploty, "C2-")
NaCsPlot.plot_survival_data(data_15m, fmt="C3s", label="Cold")
plot(fit_15m.plotx, fit_15m.ploty, "C3-")
xlim([0, 22])
ylim([0, 0.9])
legend(fontsize="small", ncol=2, borderpad=0.2, labelspacing=0.2,
       handletextpad=0.3, columnspacing=0.2, borderaxespad=0.4, loc="center left")
text(11, 0.76, "\$\\gamma_{Na}=$(fit_15na.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C0")
text(11, 0.63, "\$\\gamma_{Cs}=$(fit_15cs.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C1")
text(9, 0.32, "\$\\gamma_{Hot}=$(fit_15a.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C2")
text(9, 0.09, "\$\\gamma_{Cold}=$(fit_15m.uncs[2])\\ \\mathrm{ms}^{-1}\$",
     fontsize="small", color="C3")
text(10, 0.45, "\$\\gamma_{PA}=$(γ_15pa)\\ \\mathrm{ms}^{-1}\$",
     fontsize="small")
title("288358 GHz, 15 mW")
grid()
xlabel("Time (ms)")
ylabel("Two-body survival")

const γ_6pa = fit_6m.uncs[2] - fit_6na.uncs[2] - fit_6cs.uncs[2]
subplot(2, 2, 2)
NaCsPlot.plot_survival_data(data_6na, fmt="C0s", label="Na")
plot(fit_6na.plotx, fit_6na.ploty, "C0-")
NaCsPlot.plot_survival_data(data_6cs, fmt="C1s", label="Cs")
plot(fit_6cs.plotx, fit_6cs.ploty, "C1-")
NaCsPlot.plot_survival_data(data_6a, fmt="C2s", label="Hot")
plot(fit_6a.plotx, fit_6a.ploty, "C2-")
NaCsPlot.plot_survival_data(data_6m, fmt="C3s", label="Cold")
plot(fit_6m.plotx, fit_6m.ploty, "C3-")
xlim([0, 170])
ylim([0, 0.9])
legend(fontsize="small", ncol=2, borderpad=0.2, labelspacing=0.2,
       handletextpad=0.3, columnspacing=0.2, borderaxespad=0.4, loc="center left")
text(88, 0.76, "\$\\gamma_{Na}=$(fit_6na.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C0")
text(88, 0.63, "\$\\gamma_{Cs}=$(fit_6cs.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C1")
text(72, 0.34, "\$\\gamma_{Hot}=$(fit_6a.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C2")
text(72, 0.11, "\$\\gamma_{Cold}=$(fit_6m.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C3")
text(78, 0.45, "\$\\gamma_{PA}=$(γ_6pa * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small")
title("288358 GHz, 6 mW")
grid()
xlabel("Time (ms)")
ylabel("Two-body survival")

const γ_3pa = fit_3m.uncs[2] - fit_3cs.uncs[2] - fit_3na.uncs[2]
subplot(2, 2, 3)
NaCsPlot.plot_survival_data(data_3na, fmt="C0s", label="Na")
plot(fit_3na.plotx, fit_3na.ploty, "C0-")
NaCsPlot.plot_survival_data(data_3cs, fmt="C1s", label="Cs")
plot(fit_3cs.plotx, fit_3cs.ploty, "C1-")
NaCsPlot.plot_survival_data(data_3a, fmt="C2s", label="Hot")
plot(fit_3a.plotx, fit_3a.ploty, "C2-")
NaCsPlot.plot_survival_data(data_3m, fmt="C3s", label="Cold")
plot(fit_3m.plotx, fit_3m.ploty, "C3-")
xlim([0, 330])
ylim([0, 0.9])
legend(fontsize="small", ncol=2, borderpad=0.2, labelspacing=0.2,
       handletextpad=0.3, columnspacing=0.2, borderaxespad=0.4, loc="center left")
text(88 * 2, 0.76, "\$\\gamma_{Na}=$(fit_3na.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C0")
text(88 * 2, 0.63, "\$\\gamma_{Cs}=$(fit_3cs.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C1")
text(72 * 2, 0.34, "\$\\gamma_{Hot}=$(fit_3a.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C2")
text(72 * 2, 0.05, "\$\\gamma_{Cold}=$(fit_3m.uncs[2] * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small", color="C3")
text(160, 0.45, "\$\\gamma_{PA}=$(γ_3pa * 1000)\\ \\mathrm{s}^{-1}\$",
     fontsize="small")
title("288358 GHz, 3 mW")
grid()
xlabel("Time (ms)")
ylabel("Two-body survival")

fit_pa_1 = fit_data(gen_power_model(1.58), [3, 6, 15], [γ_3pa.a, γ_6pa.a, γ_15pa.a] .* 1000,
                    [γ_3pa.s, γ_6pa.s, γ_15pa.s] .* 1000, [1.0]; plot_lo=2.6)
fit_pa_2 = fit_data(gen_power_model(2.58), [3, 6, 15], [γ_3pa.a, γ_6pa.a, γ_15pa.a] .* 1000,
                    [γ_3pa.s, γ_6pa.s, γ_15pa.s] .* 1000, [1.0]; plot_lo=2.6)

subplot(2, 2, 4)
errorbar([3, 6, 15], [γ_3pa.a, γ_6pa.a, γ_15pa.a] .* 1000,
         [γ_3pa.s, γ_6pa.s, γ_15pa.s] .* 1000, fmt="C0s")
plot(fit_pa_1.plotx, fit_pa_1.ploty, "C1", label="\$p^{1.58}\$")
plot(fit_pa_2.plotx, fit_pa_2.ploty, "C2", label="\$p^{2.58}\$")
legend()
grid()
xscale("log")
yscale("log")
xlabel("Power (mW)")
ylabel("PA rate (\$s^{-1}\$)")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
