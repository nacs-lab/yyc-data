#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures

const expdir = joinpath(@__DIR__, "../../../experiments")

const inames = ["nacs_201808/data/data_20180815_000658.mat",
                "nacs_201808/data/data_20180815_091610.mat",
                "nacs_201808/data/data_20180815_231138.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(expdir, iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [
    OrderedDict(
        :n0=>.-(7 .+ linspace(-40, 40, 81)),
        :n1=>.-(7 .+ linspace(-40, 40, 81)),
    ),
    OrderedDict(
        :n0=>.-(7 .+ linspace(-60.0, 10.0, 81)),
        :n1=>.-(7 .+ linspace(-60.0, 10.0, 81)),
    ),
    OrderedDict(
        :cs=>.-linspace(-80.0, 80.0, 161),
        :na=>.-linspace(-80.0, 80.0, 161),
    ),
]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs_cs = select_datas(datas, NaCsData.select_single((1, 2), (4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (4,)), maxcnts, specs)
const datas_nacs_na = select_datas(datas, NaCsData.select_single((1, 2), (3,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "../figures/interaction_shift_spectrum")

data_cs = [datas_cs[1][:n0]; datas_cs[2][:n0]; datas_cs[3][:cs]]
data_nacs_cs = [datas_nacs_cs[1][:n0]; datas_nacs_cs[2][:n0]; datas_nacs_cs[3][:cs]]

data_na = datas_na[3][:na]
data_nacs_na = datas_nacs_na[3][:na]

data_nacs_cs_0 = filter(x->x < -21, data_nacs_cs)
data_nacs_cs_h = filter(x->-21 <= x < 20, data_nacs_cs)
data_nacs_cs_2 = filter(x->20 <= x < 58, data_nacs_cs)
data_nacs_cs_4 = filter(x->58 <= x, data_nacs_cs)

data_nacs_na_h = filter(x->x < 0, data_nacs_na)
data_nacs_na_0 = filter(x->0 <= x < 35, data_nacs_na)
data_nacs_na_2 = filter(x->35 <= x, data_nacs_na)

function model_gaussian(x, p)
    p[1] ./ exp.(((x .- p[2]) ./ p[3]).^2)
end

fit_cs = fit_survival(model_gaussian, filter(x->-20 < x < 25, data_cs), [0.75, 2.8, 3])
fit_nacs_cs_0 = fit_survival(model_gaussian, filter(x->-50 < x, data_nacs_cs_0),
                             [0.5, -28.4, 4])
fit_nacs_cs_h = fit_survival(model_gaussian, data_nacs_cs_h, [0.2, 3.2, 5])
fit_nacs_cs_2 = fit_survival(model_gaussian, data_nacs_cs_2, [0.15, 35, 5])
fit_na = fit_survival(model_gaussian, filter(x->-40 < x < 10, data_na), [0.75, -17, 2])
fit_nacs_na_h = fit_survival(model_gaussian, filter(x->-23 < x < -10, data_nacs_na_h),
                             [0.25, -17, 3])
fit_nacs_na_0 = fit_survival(model_gaussian, data_nacs_na_0, [0.4, 15, 3])
fit_nacs_na_2 = fit_survival(model_gaussian, filter(x->48 < x < 63, data_nacs_na_2),
                             [0.15, 35, 5])

@show fit_cs.uncs
@show fit_nacs_cs_0.uncs
@show fit_nacs_cs_h.uncs
@show fit_nacs_cs_2.uncs
@show fit_na.uncs
@show fit_nacs_na_h.uncs
@show fit_nacs_na_0.uncs
@show fit_nacs_na_2.uncs

const ptrprops = Dict(:width=>2, :headlength=>6, :headwidth=>6, :color=>"m")

figure(figsize=[6.4, 4.0])
errorbar([], [], [], fmt="C1s-", label="Cs only")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2], data_cs),
                            fmt="C1s", markersize=4)
plot(fit_cs.plotx .- fit_cs.param[2], fit_cs.ploty, "C1")
grid()
legend(fontsize="small")
text(-90, 0.792, ("    \$|\\mathrm{Cs(4,4)}\\rangle\$\n" *
                  "\$\\rightarrow|\\mathrm{Cs(3,3)}\\rangle\$"),
     ha="left", va="top", fontsize=13)
ylim([0, 0.8])
xlabel("Frequency Shift (kHz)")
ylabel("Cs Spin Flip")
NaCsPlot.maybe_save("$(prefix)_42_32_cs")

figure(figsize=[6.4, 4.0])
errorbar([], [], [], fmt="C1s-", label="Cs only")
errorbar([], [], [], fmt="C0o-", label="Na + Cs")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2], data_cs),
                            fmt="C1s", markersize=4)
plot(fit_cs.plotx .- fit_cs.param[2], fit_cs.ploty, "C1")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2],
                                                [data_nacs_cs_0; data_nacs_cs_h;
                                                 data_nacs_cs_2]),
                            fmt="C0o", markersize=4)
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2],
                                                data_nacs_cs_4),
                            fmt="C0o-", markersize=4)
plot(fit_nacs_cs_0.plotx .- fit_cs.param[2], fit_nacs_cs_0.ploty, "C0")
plot(fit_nacs_cs_h.plotx .- fit_cs.param[2], fit_nacs_cs_h.ploty, "C0")
plot(fit_nacs_cs_2.plotx .- fit_cs.param[2], fit_nacs_cs_2.ploty, "C0")
grid()
legend(fontsize="small")
text(-90, 0.792, ("    \$|\\mathrm{Na(2,2)\\ Cs(4,4)}\\rangle\$\n" *
                  "\$\\rightarrow|\\mathrm{Na(2,2)\\ Cs(3,3)}\\rangle\$"),
     ha="left", va="top", fontsize=13)
ylim([0, 0.8])
xlabel("Frequency Shift (kHz)")
ylabel("Cs Spin Flip")
NaCsPlot.maybe_save("$(prefix)_42_32_nostate")

figure(figsize=[6.4, 4.0])
errorbar([], [], [], fmt="C1s-", label="Cs only")
errorbar([], [], [], fmt="C0o-", label="Na + Cs")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2], data_cs),
                            fmt="C1s", markersize=4)
plot(fit_cs.plotx .- fit_cs.param[2], fit_cs.ploty, "C1")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2],
                                                [data_nacs_cs_0; data_nacs_cs_h;
                                                 data_nacs_cs_2]),
                            fmt="C0o", markersize=4)
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs.param[2],
                                                data_nacs_cs_4),
                            fmt="C0o-", markersize=4)
plot(fit_nacs_cs_0.plotx .- fit_cs.param[2], fit_nacs_cs_0.ploty, "C0")
plot(fit_nacs_cs_h.plotx .- fit_cs.param[2], fit_nacs_cs_h.ploty, "C0")
plot(fit_nacs_cs_2.plotx .- fit_cs.param[2], fit_nacs_cs_2.ploty, "C0")
grid()
legend(fontsize="small")
text(-90, 0.792, ("    \$|\\mathrm{Na(2,2)\\ Cs(4,4)}\\rangle\$\n" *
                  "\$\\rightarrow|\\mathrm{Na(2,2)\\ Cs(3,3)}\\rangle\$"),
     ha="left", va="top", fontsize=13)
x_cs_0 = fit_nacs_cs_0.param[2] - fit_cs.param[2]
annotate("\$|0,0\\rangle\$", xy=(x_cs_0, 0.54),
         xytext=(x_cs_0, 0.60), arrowprops=ptrprops, fontsize=14, ha="center", color="m")
x_cs_2 = fit_nacs_cs_2.param[2] - fit_cs.param[2]
annotate("\$|2,0\\rangle\$", xy=(x_cs_2, 0.20),
         xytext=(x_cs_2, 0.26), arrowprops=ptrprops, fontsize=14, ha="center", color="m")
x_cs_4 = 76.2
annotate("\$|4,0\\rangle\$", xy=(x_cs_4, 0.07),
         xytext=(x_cs_4, 0.13), arrowprops=ptrprops, fontsize=14, ha="center", color="m")
ylim([0, 0.8])
xlabel("Frequency Shift (kHz)")
ylabel("Cs Spin Flip")
NaCsPlot.maybe_save("$(prefix)_42_32")

figure(figsize=[6.4, 4.0])
errorbar([], [], [], fmt="C1s-", label="Na only")
errorbar([], [], [], fmt="C0o-", label="Na + Cs")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_na.param[2], data_na),
                            fmt="C1s", markersize=4)
plot(fit_na.plotx .- fit_na.param[2], fit_na.ploty, "C1")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_na.param[2], data_nacs_na),
                            fmt="C0o", markersize=4)
plot(fit_nacs_na_h.plotx .- fit_na.param[2], fit_nacs_na_h.ploty, "C0")
plot(fit_nacs_na_0.plotx .- fit_na.param[2], fit_nacs_na_0.ploty, "C0")
plot(fit_nacs_na_2.plotx .- fit_na.param[2], fit_nacs_na_2.ploty, "C0")
grid()
legend(fontsize="small")
text(-70, 0.792, ("    \$|\\mathrm{Na(2,2)\\ Cs(3,3)}\\rangle\$\n" *
                  "\$\\rightarrow|\\mathrm{Na(1,1)\\ Cs(3,3)}\\rangle\$"),
     ha="left", va="top", fontsize=13)
x_na_0 = fit_nacs_na_0.param[2] - fit_na.param[2]
annotate("\$|0,0\\rangle\$", xy=(x_na_0, 0.45),
         xytext=(x_na_0, 0.51), arrowprops=ptrprops, fontsize=14, ha="center", color="m")
x_na_2 = fit_nacs_na_2.param[2] - fit_na.param[2]
annotate("\$m_{\\mathrm{rel},z}\\!+\\!m_{\\mathrm{COM},z}\\!=\\!2\$", xy=(x_na_2, 0.11),
         xytext=(x_na_2, 0.17), arrowprops=ptrprops, fontsize=14, ha="center", color="m")
ylim([0, 0.8])
xlabel("Frequency Shift (kHz)")
ylabel("Na Spin Flip")
NaCsPlot.maybe_save("$(prefix)_32_31")

NaCsPlot.maybe_show()
