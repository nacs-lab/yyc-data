#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200109_121138.mat",
                "data_20200109_150852.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [368.1571 .+ (-30:5:30) .* 1e-3,
               367.5 .+ (-40:5:40) .* 1e-3]
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
const data_460 = datas_nacs[1]
const data_550 = datas_nacs[2]

const prefix = joinpath(@__DIR__, "imgs", "data_20200109_raman_det_red_vs_blue_3311")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
const fit_460 = fit_survival(model_lorentzian, data_460, [0.7, 0.3, 368.157, 0.05])
const fit_550 = fit_survival(model_lorentzian, data_550, [0.7, 0.3, 367.5, 0.05])

const data_460′ = NaCsData.map_params((i, v) -> (v .- fit_460.param[3]) .* 1000, data_460)
const data_550′ = NaCsData.map_params((i, v) -> (v .- fit_550.param[3]) .* 1000, data_550)

figure()
NaCsPlot.plot_survival_data(data_460′, fmt="C1.", label="306460.7 GHz")
plot((fit_460.plotx .- fit_460.param[3]) .* 1000, fit_460.ploty, "C1-")
NaCsPlot.plot_survival_data(data_550′, fmt="C0.", label="306550.7 GHz")
plot((fit_550.plotx .- fit_550.param[3]) .* 1000, fit_550.ploty, "C0-")
title("8 mW, 5 ms")
# # ylim([0.18, 0.42])
text(-52, 0.51, "\$f=$(fit_460.uncs[3])\$ MHz", color="C1", fontsize="small")
text(-52, 0.54, "\$\\Gamma=$(fit_460.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
text(-7, 0.45, "\$f=$(fit_550.uncs[3])\$ MHz", color="C0", fontsize="small")
text(2, 0.48, "\$\\Gamma=$(fit_550.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
legend(fontsize="small", loc="upper center", ncol=2)
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
