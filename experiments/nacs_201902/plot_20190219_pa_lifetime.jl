#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190219_151703.mat",
                "data_20190220_094017.mat",
                "data_20190220_144236.mat",
                "data_20190220_153441.mat",
                "data_20190220_171238.mat",
                "data_20190220_215125.mat", # 682 GHz
                "data_20190221_081927.mat",
                "data_20190221_143511.mat",
                "data_20190221_154307.mat",
                "data_20190221_223037.mat",
                "data_20190222_110328.mat", # 691.5 GHz
                "data_20190222_185937.mat",
                "data_20190223_001919.mat",
                "data_20190223_090035.mat",
                "data_20190223_123333.mat",
                "data_20190223_140914.mat",
                "data_20190224_084424.mat",
                "data_20190224_163643.mat", # 697.5 GHz
                "data_20190224_213643.mat",
                "data_20190225_065431.mat",
                "data_20190225_101432.mat",
                "data_20190225_133508.mat",
                "data_20190225_185048.mat",
                "data_20190225_222710.mat", # 705 GHz
                "data_20190226_020710.mat", # 707 GHz
                "data_20190226_054558.mat", # 709 GHz
                "data_20190226_094328.mat", # 709 GHz
                "data_20190226_131609.mat", # 707 GHz
                "data_20190226_155840.mat",
                "data_20190226_205050.mat",
                "data_20190227_013017.mat",
                "data_20190227_092022.mat",
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int), # 682 GHz
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int), # 691.5 GHz
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int), # 697.5
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int), # 705 GHz
                 typemax(Int), # 707 GHz
                 typemax(Int), # 709 GHz
                 typemax(Int), # 709 GHz
                 typemax(Int), # 707 GHz
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 ]
const specs = [([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500], # 682 GHz
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 5, 10, 20, 50, 100, 200], # 691.5 GHz
                [0, 5, 10, 20, 50, 100, 200]),
               [0, 1, 2, 4, 6, 8, 10, 12, 15, 20],
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               [0, 2, 5, 10, 20, 50, 100],
               [0, 1, 2, 3, 4, 6, 8],
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100], # 697.5
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]),
               ([0, 5, 10, 20, 50, 100, 200], # 705 GHz
                [0, 1, 2, 4, 6, 8, 10, 15, 20, 30]),
               ([0, 5, 10, 20, 50, 100, 200], # 707 GHz
                [0, 1, 2, 4, 6, 8, 10, 15, 20, 30]),
               ([0, 5, 10, 20, 50, 100, 200], # 709 GHz
                [0, 1, 2, 4, 6, 8, 10, 15, 20, 30]),
               ([500], # 709 GHz
                [0, 20, 30, 50, 100]),
               [0, 20, 30, 50, 100], # 709 GHz
               ([0, 5, 10, 20, 50, 100, 200, 500], # 712 GHz
                [0, 5, 10, 20, 50, 100, 200, 500]),
               ([0, 5, 10, 20, 50, 100, 200, 500], # 715 GHz
                [0, 5, 10, 20, 50, 100, 200, 500]),
               ([0, 5, 10, 20, 50, 100, 200, 500], # 718 GHz
                [0, 5, 10, 20, 50, 100, 200, 500]),
               ([0, 5, 10, 20, 50, 100, 200, 500], # 722 GHz
                [0, 5, 10, 20, 50, 100, 200, 500]),
                ]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function fit_survival(model, data, p0; plotx=nothing, plot_lo=nothing, plot_hi=nothing,
                      use_unc=true, plot_scale=1.1)
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
        fit = curve_fit(model, params, ratios[:, 2], uncs[:, 2].^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4)), maxcnts, specs)

model_exp(x, p) = p[1] .* exp.(x ./ -p[2])
model_exp2(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3] .* exp.(x ./ -p[4])

data_640 = datas_nacs[1]
data_656 = datas_nacs[2]
data_670 = [[datas_nacs[3][1]; datas_nacs[4][1]; datas_nacs[5][1]],
            [datas_nacs[3][2]; datas_nacs[4][2]; datas_nacs[5][2]]]
data_682 = datas_nacs[6]
data_688 = [[datas_nacs[7][1]; datas_nacs[8][1]; datas_nacs[9][1]],
            [datas_nacs[7][2]; datas_nacs[8][2]; datas_nacs[9][2]]]
data_690 = datas_nacs[10]
data_691_5 = [datas_nacs[11][1], [datas_nacs[11][2]; datas_nacs[12]]]
data_692_5 = datas_nacs[13]
data_693_5 = datas_nacs[14]
data_694_5 = datas_nacs[15]
data_695_5 = datas_nacs[16]
data_696_5 = datas_nacs[17]
data_697_5 = datas_nacs[18]
data_698_5 = datas_nacs[19]
data_699_5 = datas_nacs[20]
data_700_5 = datas_nacs[21]
data_701_5 = datas_nacs[22]
data_703 = datas_nacs[23]
data_705 = datas_nacs[24]
data_707 = [datas_nacs[25][1], [datas_nacs[25][2]; datas_nacs[28]]]
data_709 = [[datas_nacs[26][1]; datas_nacs[27][1]],
            [datas_nacs[26][2]; datas_nacs[27][2]]]
data_712 = datas_nacs[29]
data_715 = datas_nacs[30]
data_718 = datas_nacs[31]
data_722 = datas_nacs[32]

fit_640_hot = fit_survival(model_exp, data_640[2], [0.75, 800])
fit_640_cold = fit_survival(model_exp2, data_640[4], [0.5, 200, 0.2, 800])
fit_656_hot = fit_survival(model_exp, data_656[1], [0.75, 800])
fit_656_cold = fit_survival(model_exp2, data_656[2], [0.5, 200, 0.2, 800])
fit_670_hot = fit_survival(model_exp, data_670[1], [0.75, 800])
fit_670_cold = fit_survival(model_exp2, data_670[2], [0.5, 200, 0.2, 800])
fit_682_hot = fit_survival(model_exp, data_682[1], [0.75, 800])
fit_682_cold = fit_survival(model_exp2, data_682[2], [0.5, 200, 0.2, 800])
fit_688_hot = fit_survival(model_exp, data_688[1], [0.75, 800])
fit_688_cold = fit_survival(model_exp2, data_688[2], [0.5, 200, 0.2, 800])
fit_690_hot = fit_survival(model_exp, data_690[1], [0.75, 800])
fit_690_cold = fit_survival(model_exp2, data_690[2], [0.5, 200, 0.2, 800])
fit_691_5_hot = fit_survival(model_exp, data_691_5[1], [0.75, 800])
fit_691_5_cold = fit_survival(model_exp2, data_691_5[2], [0.5, 200, 0.2, 800])
fit_692_5_hot = fit_survival(model_exp, data_692_5[1], [0.75, 70])
fit_692_5_cold = fit_survival(model_exp, data_692_5[2], [0.2, 40])
fit_693_5_hot = fit_survival(model_exp, data_693_5, [0.75, 60])
fit_694_5_hot = fit_survival(model_exp, data_694_5, [0.75, 5])
fit_695_5_hot = fit_survival(model_exp, data_695_5[1], [0.75, 100])
fit_695_5_cold = fit_survival(model_exp, data_695_5[2], [0.2, 40])
fit_696_5_hot = fit_survival(model_exp, data_696_5[1], [0.75, 100])
fit_696_5_cold = fit_survival(model_exp, data_696_5[2], [0.2, 40])
fit_697_5_hot = fit_survival(model_exp, data_697_5[1], [0.75, 100])
fit_697_5_cold = fit_survival(model_exp, data_697_5[2], [0.2, 40])
fit_698_5_hot = fit_survival(model_exp, data_698_5[1], [0.75, 100])
fit_698_5_cold = fit_survival(model_exp, data_698_5[2], [0.2, 40])
fit_699_5_hot = fit_survival(model_exp, data_699_5[1], [0.75, 100])
fit_699_5_cold = fit_survival(model_exp, data_699_5[2], [0.2, 40])
fit_700_5_hot = fit_survival(model_exp, data_700_5[1], [0.75, 100])
fit_700_5_cold = fit_survival(model_exp, data_700_5[2], [0.2, 40])
fit_701_5_hot = fit_survival(model_exp, data_701_5[1], [0.75, 100])
fit_701_5_cold = fit_survival(model_exp, data_701_5[2], [0.2, 40])
fit_703_hot = fit_survival(model_exp, data_703[1], [0.75, 100])
fit_703_cold = fit_survival(model_exp2, data_703[2], [0.3, 10, 0.1, 100])
fit_705_hot = fit_survival(model_exp, data_705[1], [0.75, 100])
fit_705_cold = fit_survival(model_exp2, data_705[2], [0.3, 10, 0.1, 100])
fit_707_hot = fit_survival(model_exp, data_707[1], [0.75, 200])
fit_707_cold = fit_survival(model_exp2, data_707[2], [0.5, 40, 0.2, 200])
fit_709_hot = fit_survival(model_exp, data_709[1], [0.75, 200])
fit_709_cold = fit_survival(model_exp2, data_709[2], [0.5, 40, 0.2, 200])
fit_712_hot = fit_survival(model_exp, data_712[1], [0.75, 500])
fit_712_cold = fit_survival(model_exp2, data_712[2], [0.5, 100, 0.2, 500])
fit_715_hot = fit_survival(model_exp, data_715[1], [0.75, 500])
fit_715_cold = fit_survival(model_exp2, data_715[2], [0.5, 100, 0.2, 500])
fit_718_hot = fit_survival(model_exp, data_718[1], [0.75, 500])
fit_718_cold = fit_survival(model_exp2, data_718[2], [0.5, 100, 0.2, 500])
fit_722_hot = fit_survival(model_exp, data_722[1], [0.75, 500])
fit_722_cold = fit_survival(model_exp2, data_722[2], [0.5, 100, 0.2, 500])

struct Line
    freqs::Vector{Float64}
    p::Vector{Float64}
    p_s::Vector{Float64}
    r::Vector{Float64}
    r_s::Vector{Float64}
    Line() = new(Float64[], Float64[], Float64[], Float64[], Float64[])
end

const line_hot = Line()
const line_cold1 = Line()
const line_cold2 = Line()

function add_point!(line::Line, freq, fit, i0=1)
    push!(line.freqs, freq)
    push!(line.p, fit.param[i0])
    push!(line.p_s, fit.unc[i0])
    push!(line.r, 1000 / fit.param[i0 + 1])
    push!(line.r_s, fit.unc[i0 + 1] / fit.param[i0 + 1]^2 * 1000)
end

function add_fit!(freq, hot, cold=nothing)
    add_point!(line_hot, freq, hot)
    cold === nothing && return
    if length(cold.param) == 2
        add_point!(line_cold1, freq, cold)
    elseif cold.param[4] < cold.param[2]
        add_point!(line_cold1, freq, cold, 1)
        add_point!(line_cold2, freq, cold, 3)
    else
        add_point!(line_cold2, freq, cold, 1)
        add_point!(line_cold1, freq, cold, 3)
    end
end

add_fit!(640, fit_640_hot, fit_640_cold)
add_fit!(656, fit_656_hot, fit_656_cold)
add_fit!(670, fit_670_hot, fit_670_cold)
add_fit!(682, fit_682_hot, fit_682_cold)
add_fit!(688, fit_688_hot, fit_688_cold)
add_fit!(690, fit_690_hot, fit_690_cold)
add_fit!(691.5, fit_691_5_hot, fit_691_5_cold)
add_fit!(692.5, fit_692_5_hot, fit_692_5_cold)
add_fit!(693.5, fit_693_5_hot)
add_fit!(694.5, fit_694_5_hot)
add_fit!(695.5, fit_695_5_hot, fit_695_5_cold)
add_fit!(696.5, fit_696_5_hot, fit_696_5_cold)
add_fit!(697.5, fit_697_5_hot, fit_697_5_cold)
add_fit!(698.5, fit_698_5_hot, fit_698_5_cold)
add_fit!(699.5, fit_699_5_hot, fit_699_5_cold)
add_fit!(700.5, fit_700_5_hot, fit_700_5_cold)
add_fit!(701.5, fit_701_5_hot, fit_701_5_cold)
add_fit!(703, fit_703_hot, fit_703_cold)
add_fit!(705, fit_705_hot, fit_705_cold)
add_fit!(707, fit_707_hot, fit_707_cold)
add_fit!(709, fit_709_hot, fit_709_cold)
add_fit!(712, fit_712_hot, fit_712_cold)
add_fit!(715, fit_715_hot, fit_715_cold)
add_fit!(718, fit_718_hot, fit_718_cold)
add_fit!(722, fit_722_hot, fit_722_cold)

const prefix = joinpath(@__DIR__, "imgs", "data_20190219_ps_lifetime")

figure()
NaCsPlot.plot_survival_data(data_640[2], fmt="C0.", label="640")
plot(fit_640_hot.plotx, fit_640_hot.ploty, "C0")
NaCsPlot.plot_survival_data(data_656[1], fmt="C1.", label="656")
plot(fit_656_hot.plotx, fit_656_hot.ploty, "C1")
NaCsPlot.plot_survival_data(data_670[1], fmt="C2.", label="670")
plot(fit_670_hot.plotx, fit_670_hot.ploty, "C2")
NaCsPlot.plot_survival_data(data_682[1], fmt="C3.", label="682")
plot(fit_682_hot.plotx, fit_682_hot.ploty, "C3")
NaCsPlot.plot_survival_data(data_688[1], fmt="C4.", label="688")
plot(fit_688_hot.plotx, fit_688_hot.ploty, "C4")
NaCsPlot.plot_survival_data(data_690[1], fmt="C5.", label="690")
plot(fit_690_hot.plotx, fit_690_hot.ploty, "C5")
NaCsPlot.plot_survival_data(data_691_5[1], fmt="C6.", label="691.5")
plot(fit_691_5_hot.plotx, fit_691_5_hot.ploty, "C6")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 550])
title("2-Body Lifetime (Hot)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_hot1")

figure()
NaCsPlot.plot_survival_data(data_640[4], fmt="C0.", label="640")
plot(fit_640_cold.plotx, fit_640_cold.ploty, "C0")
NaCsPlot.plot_survival_data(data_656[2], fmt="C1.", label="656")
plot(fit_656_cold.plotx, fit_656_cold.ploty, "C1")
NaCsPlot.plot_survival_data(data_670[2], fmt="C2.", label="670")
plot(fit_670_cold.plotx, fit_670_cold.ploty, "C2")
NaCsPlot.plot_survival_data(data_682[2], fmt="C3.", label="682")
plot(fit_682_cold.plotx, fit_682_cold.ploty, "C3")
NaCsPlot.plot_survival_data(data_688[2], fmt="C4.", label="688")
plot(fit_688_cold.plotx, fit_688_cold.ploty, "C4")
NaCsPlot.plot_survival_data(data_690[2], fmt="C5.", label="690")
plot(fit_690_cold.plotx, fit_690_cold.ploty, "C5")
NaCsPlot.plot_survival_data(data_691_5[2], fmt="C6.", label="691.5")
plot(fit_691_5_cold.plotx, fit_691_5_cold.ploty, "C6")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 250])
title("2-Body Lifetime (Cold)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cold1")

figure()
NaCsPlot.plot_survival_data(data_692_5[1], fmt="C0.", label="692.5")
plot(fit_692_5_hot.plotx, fit_692_5_hot.ploty, "C0")
NaCsPlot.plot_survival_data(data_693_5, fmt="C1.", label="693.5")
plot(fit_693_5_hot.plotx, fit_693_5_hot.ploty, "C1")
NaCsPlot.plot_survival_data(data_694_5, fmt="C2.", label="694.5")
plot(fit_694_5_hot.plotx, fit_694_5_hot.ploty, "C2")
NaCsPlot.plot_survival_data(data_695_5[1], fmt="C3.", label="695.5")
plot(fit_695_5_hot.plotx, fit_695_5_hot.ploty, "C3")
NaCsPlot.plot_survival_data(data_696_5[1], fmt="C4.", label="696.5")
plot(fit_696_5_hot.plotx, fit_696_5_hot.ploty, "C4")
NaCsPlot.plot_survival_data(data_697_5[1], fmt="C5.", label="697.5")
plot(fit_697_5_hot.plotx, fit_697_5_hot.ploty, "C5")
NaCsPlot.plot_survival_data(data_698_5[1], fmt="C6.", label="698.5")
plot(fit_698_5_hot.plotx, fit_698_5_hot.ploty, "C6")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 105])
title("2-Body Lifetime (Hot)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_hot2")

figure()
NaCsPlot.plot_survival_data(data_692_5[2], fmt="C0.", label="692.5")
plot(fit_692_5_cold.plotx, fit_692_5_cold.ploty, "C0")
NaCsPlot.plot_survival_data(data_695_5[2], fmt="C3.", label="695.5")
plot(fit_695_5_cold.plotx, fit_695_5_cold.ploty, "C3")
NaCsPlot.plot_survival_data(data_696_5[2], fmt="C4.", label="696.5")
plot(fit_696_5_cold.plotx, fit_696_5_cold.ploty, "C4")
NaCsPlot.plot_survival_data(data_697_5[2], fmt="C5.", label="697.5")
plot(fit_697_5_cold.plotx, fit_697_5_cold.ploty, "C5")
NaCsPlot.plot_survival_data(data_698_5[2], fmt="C6.", label="698.5")
plot(fit_698_5_cold.plotx, fit_698_5_cold.ploty, "C6")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 0.25])
xlim([0, 22])
title("2-Body Lifetime (Cold)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cold2")

figure()
NaCsPlot.plot_survival_data(data_699_5[1], fmt="C0.", label="699.5")
plot(fit_699_5_hot.plotx, fit_699_5_hot.ploty, "C0")
NaCsPlot.plot_survival_data(data_700_5[1], fmt="C1.", label="700.5")
plot(fit_700_5_hot.plotx, fit_700_5_hot.ploty, "C1")
NaCsPlot.plot_survival_data(data_701_5[1], fmt="C2.", label="701.5")
plot(fit_701_5_hot.plotx, fit_701_5_hot.ploty, "C2")
NaCsPlot.plot_survival_data(data_703[1], fmt="C3.", label="703")
plot(fit_703_hot.plotx, fit_703_hot.ploty, "C3")
NaCsPlot.plot_survival_data(data_705[1], fmt="C4.", label="705")
plot(fit_705_hot.plotx, fit_705_hot.ploty, "C4")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 105])
title("2-Body Lifetime (Hot)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_hot3")

figure()
NaCsPlot.plot_survival_data(data_699_5[2], fmt="C0.", label="699.5")
plot(fit_699_5_cold.plotx, fit_699_5_cold.ploty, "C0")
NaCsPlot.plot_survival_data(data_700_5[2], fmt="C1.", label="700.5")
plot(fit_700_5_cold.plotx, fit_700_5_cold.ploty, "C1")
NaCsPlot.plot_survival_data(data_701_5[2], fmt="C2.", label="701.5")
plot(fit_701_5_cold.plotx, fit_701_5_cold.ploty, "C2")
NaCsPlot.plot_survival_data(data_703[2], fmt="C3.", label="703")
plot(fit_703_cold.plotx, fit_703_cold.ploty, "C3")
NaCsPlot.plot_survival_data(data_705[2], fmt="C4.", label="705")
plot(fit_705_cold.plotx, fit_705_cold.ploty, "C4")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 0.6])
xlim([0, 22])
title("2-Body Lifetime (Cold)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cold3")

figure()
NaCsPlot.plot_survival_data(data_707[1], fmt="C0.", label="707")
plot(fit_707_hot.plotx, fit_707_hot.ploty, "C0")
NaCsPlot.plot_survival_data(data_709[1], fmt="C1.", label="709")
plot(fit_709_hot.plotx, fit_709_hot.ploty, "C1")
NaCsPlot.plot_survival_data(data_712[1], fmt="C2.", label="712")
plot(fit_712_hot.plotx, fit_712_hot.ploty, "C2")
NaCsPlot.plot_survival_data(data_715[1], fmt="C3.", label="715")
plot(fit_715_hot.plotx, fit_715_hot.ploty, "C3")
NaCsPlot.plot_survival_data(data_718[1], fmt="C4.", label="718")
plot(fit_718_hot.plotx, fit_718_hot.ploty, "C4")
NaCsPlot.plot_survival_data(data_722[1], fmt="C5.", label="722")
plot(fit_722_hot.plotx, fit_722_hot.ploty, "C5")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 550])
title("2-Body Lifetime (Hot)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_hot4")

figure()
NaCsPlot.plot_survival_data(data_707[2], fmt="C0.", label="707")
plot(fit_707_cold.plotx, fit_707_cold.ploty, "C0")
NaCsPlot.plot_survival_data(data_709[2], fmt="C1.", label="709")
plot(fit_709_cold.plotx, fit_709_cold.ploty, "C1")
NaCsPlot.plot_survival_data(data_712[2], fmt="C2.", label="712")
plot(fit_712_cold.plotx, fit_712_cold.ploty, "C2")
NaCsPlot.plot_survival_data(data_715[2], fmt="C3.", label="715")
plot(fit_715_cold.plotx, fit_715_cold.ploty, "C3")
NaCsPlot.plot_survival_data(data_718[2], fmt="C4.", label="718")
plot(fit_718_cold.plotx, fit_718_cold.ploty, "C4")
NaCsPlot.plot_survival_data(data_722[2], fmt="C5.", label="722")
plot(fit_722_cold.plotx, fit_722_cold.ploty, "C5")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 250])
title("2-Body Lifetime (Cold)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cold4")

figure()
errorbar(line_hot.freqs, line_hot.p, line_hot.p_s, fmt=".-", label="Hot")
errorbar(line_cold1.freqs, line_cold1.p, line_cold1.p_s, fmt=".-", label="Cold (slow)")
errorbar(line_cold2.freqs, line_cold2.p, line_cold2.p_s, fmt=".-", label="Cold (fast)")
grid()
legend()
ylim([0, 1])
title("2-Body lifetime fits")
xlabel("Frequency (288XXX GHz)")
ylabel("Initial survival")
NaCsPlot.maybe_save("$(prefix)_p0s")

figure()
errorbar(line_hot.freqs, line_hot.r, line_hot.r_s, fmt=".-", label="Hot")
errorbar(line_cold1.freqs, line_cold1.r, line_cold1.r_s, fmt=".-", label="Cold (slow)")
errorbar(line_cold2.freqs, line_cold2.r, line_cold2.r_s, fmt=".-", label="Cold (fast)")
grid()
legend()
ax = gca()
ax[:set_yscale]("log", nonposy="clip")
ax[:set_yticks]([1, 10, 100])
ax[:set_yticklabels](["1", "10", "100"])
ylim([0.7, 700])
title("2-Body lifetime fits")
xlabel("Frequency (288XXX GHz)")
ylabel("Decay rate (\$s^{-1}\$)")
NaCsPlot.maybe_save("$(prefix)_rates")

NaCsPlot.maybe_show()
