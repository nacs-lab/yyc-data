#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200108_233851.mat",
                "data_20200109_002814.mat",
                "data_20200109_093029.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [368.433 .+ (-40:5:40) .* 1e-3,
               367.22 .+ (-120:5:120) .* 1e-3,
               368.433 .+ (-40:5:40) .* 1e-3]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

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

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_480 = [datas_nacs[1]; datas_nacs[3]]
const data_530 = datas_nacs[2]

const prefix = joinpath(@__DIR__, "imgs", "data_20200108_raman_det_red_vs_blue_3311")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
const fit_480 = fit_survival(model_lorentzian, data_480, [0.7, 0.3, 368.433, 0.05])
const fit_530 = fit_survival(model_lorentzian, data_530, [0.7, 0.3, 367.22, 0.05])

const data_480′ = NaCsData.map_params((i, v) -> (v .- fit_480.param[3]) .* 1000, data_480)
const data_530′ = NaCsData.map_params((i, v) -> (v .- fit_530.param[3]) .* 1000, data_530)

figure()
NaCsPlot.plot_survival_data(data_480′, fmt="C1.", label="306480.7 GHz")
plot((fit_480.plotx .- fit_480.param[3]) .* 1000, fit_480.ploty, "C1-")
NaCsPlot.plot_survival_data(data_530′, fmt="C0.", label="306530.7 GHz")
plot((fit_530.plotx .- fit_530.param[3]) .* 1000, fit_530.ploty, "C0-")
title("8 mW, 5 ms")
# # ylim([0.18, 0.42])
text(-136, 0.3, "\$f=$(fit_480.uncs[3])\$ MHz", color="C1", fontsize="small")
text(-136, 0.335, "\$\\Gamma=$(fit_480.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
text(16, 0.47, "\$f=$(fit_530.uncs[3])\$ MHz", color="C0", fontsize="small")
text(28, 0.505, "\$\\Gamma=$(fit_530.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
legend(fontsize="small", loc="upper center", ncol=2)
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
