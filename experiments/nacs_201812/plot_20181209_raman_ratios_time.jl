#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181209_141239.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [OrderedDict(:tr50=>[0.1, 0.5, 1, 1.5, 2, 4, 7, 10, 15, 25, 40],
                           :tr200=>[0.1, 0.5, 1, 1.5, 2, 4, 7, 10, 15, 25, 40] .* 4,
                           :r50=>302.125 .+ [-300; -250; -200; -150; -120;
                                             -100:10:100;
                                             120; 150; 200; 250; 300] .* 1e-3,
                           :r200=>302.125 .+ [-300; -250; -200; -150; -120;
                                              -100:10:100;
                                              120; 150; 200; 250; 300] .* 1e-3)]
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
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
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

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2,), (3, 4,)), maxcnts, specs)

data_tr50 = datas_nacs[1][:tr50]
data_tr200 = datas_nacs[1][:tr200]
data_r50 = datas_nacs[1][:r50]
data_r200 = datas_nacs[1][:r200]

function gen_model(v0, v1)
    dv = v1 - v0
    return (x, p)->v0 .+ dv .* exp.(-p[1] ./ (1 .+ 4 .* (x .- p[2]).^2 ./ p[3]^2))
end

function model′(x, p)
    return p[4] .- p[1] ./ (1 .+ 4 .* (x .- p[2]).^2 ./ p[3]^2)
end
model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]
model_exp2(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3] .* exp.(x ./ -p[4])

fit_r50 = fit_survival(gen_model(0.25, 0.78), data_r50, [0.7, 302.125, 0.1])
@show fit_r50.uncs
fit_r200 = fit_survival(gen_model(0.25, 0.78), data_r200, [0.7, 302.125, 0.1])
@show fit_r200.uncs

fit_r50′ = fit_survival(model′, data_r50, [0.2, 302.125, 0.1, 0.7])
@show fit_r50′.uncs
fit_r200′ = fit_survival(model′, data_r200, [0.2, 302.125, 0.1, 0.7])
@show fit_r200′.uncs

fit_tr50 = fit_survival(model_exp2, data_tr50, [0.53, 1, 0.22, 10], plot_lo=0)
fit_tr200 = fit_survival(model_exp2, data_tr200, [0.5, 4, 0.22, 40], plot_lo=0)

const prefix = joinpath(@__DIR__, "imgs", "data_20181209_raman_ratios_time")

figure()
NaCsPlot.plot_survival_data(data_r50, fmt="C0.", label="0.02")
plot(fit_r50.plotx, fit_r50.ploty, "C0-")
plot(fit_r50′.plotx, fit_r50′.ploty, "C0--")
NaCsPlot.plot_survival_data(data_r200, fmt="C1.", label="0.005")
plot(fit_r200.plotx, fit_r200.ploty, "C1-")
plot(fit_r200′.plotx, fit_r200′.ploty, "C1--")
grid()
legend()
xlabel("Raman frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_freq")

figure()
NaCsPlot.plot_survival_data(data_tr50, fmt="C0.", label="0.02")
plot(fit_tr50.plotx, fit_tr50.ploty, "C0-")
NaCsPlot.plot_survival_data(data_tr200, fmt="C1.", label="0.005")
plot(fit_tr200.plotx, fit_tr200.ploty, "C1-")
grid()
ylim([0, 0.8])
xlim([0, 170])
text(8, 0.02, "\$\\tau=$(fit_tr50.uncs[2])\$ms", color="C0")
text(40, 0.2, "\$\\tau=$(fit_tr200.uncs[2])\$ms", color="C1")
legend(fontsize="small", labelspacing=0.2, borderpad=0.2,
       handletextpad=0.2, columnspacing=0.3, borderaxespad=0.2)
xlabel("Raman time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time")

NaCsPlot.maybe_show()
