#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200206_173305.mat",
                "data_20200206_204557.mat",
                "data_20200206_190515.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [linspace(0, 0.05, 21),
               0.5 .+ linspace(0, 0.05, 21),
               1 .+ linspace(0, 0.05, 21)]
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

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200206_coprop_raman_coherence")

model_sin(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4])
model_sinexp(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4]) .* exp.(.-x ./ p[5])

fit_na = fit_survival(model_sin, datas_na[1], [0.4, 0.4, 25, 0.5])
fit_na = fit_survival(model_sinexp, [datas_na[1]; datas_na[2]], [fit_na.param; 0.5])
fit_na = fit_survival(model_sinexp, [datas_na[1]; datas_na[2]; datas_na[3]], fit_na.param)

fit_cs = fit_survival(model_sin, datas_cs[1], [0.4, 0.4, 30, 0.5])
fit_cs = fit_survival(model_sinexp, [datas_cs[1]; datas_cs[2]], [fit_cs.param; 0.5])
fit_cs = fit_survival(model_sinexp, [datas_cs[1]; datas_cs[2]; datas_cs[3]], fit_cs.param)

# @show fit.uncs

f, axs = subplots(1, 3, sharey=true, facecolor="w", figsize=(6.5, 4.8))
subplots_adjust(wspace=0.05, hspace=0.05)
d = 0.015 # how big to make the diagonal lines in axes coordinates

sca(axs[1])
axs[1][:spines]["right"][:set_visible](false)
axs[1][:yaxis][:tick_left]()
axs[1][:tick_params](labelright=false)
NaCsPlot.plot_survival_data(datas_na[1], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_cs[1], fmt="C1.")
plot(fit_cs.plotx, fit_cs.ploty, "C1-")
grid()
xlim([0, 0.06])
ylim([0, 1])
kwargs = Dict(:transform=>axs[1][:transAxes], :color=>"k", :clip_on=>false)
axs[1][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[1][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)
ylabel("Na Survival")

sca(axs[2])
axs[2][:spines]["left"][:set_visible](false)
axs[2][:spines]["right"][:set_visible](false)
axs[2][:yaxis][:set_ticks_position]("none")
axs[2][:tick_params](labelright=false)
NaCsPlot.plot_survival_data(datas_na[2], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_cs[2], fmt="C1.")
plot(fit_cs.plotx, fit_cs.ploty, "C1-")
grid()
xlim([0.495, 0.555])
ylim([0, 1])
kwargs = Dict(:transform=>axs[2][:transAxes], :color=>"k", :clip_on=>false)
axs[2][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[2][:plot]((-d, d), (-d, d); kwargs...)
axs[2][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[2][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[3])
axs[3][:spines]["left"][:set_visible](false)
axs[3][:yaxis][:tick_right]()
NaCsPlot.plot_survival_data(datas_na[3], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_cs[3], fmt="C1.")
plot(fit_cs.plotx, fit_cs.ploty, "C1-")
grid()
xlim([0.995, 1.055])
ylim([0, 1])
kwargs = Dict(:transform=>axs[3][:transAxes], :color=>"k", :clip_on=>false)
axs[3][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[3][:plot]((-d, d), (-d, d); kwargs...)

f[:add_subplot](111, frameon=false)
tick_params(labelcolor="none", top=false, bottom=false, left=false, right=false)
xlabel("Raman time (ms)")
xlim([0, 1])
ylim([0, 1])

text(0.42, 0.02, "\$\\nu_{Na}=$(fit_na.uncs[3])\$ kHz", color="C0")
text(0.50, 0.105, "\$\\tau_{Na}=$(fit_na.uncs[5])\$ ms", color="C0")
text(0.41, 0.72, "\$\\nu_{Cs}=$(fit_cs.uncs[3])\$ kHz", color="C1")
text(0.46, 0.635, "\$\\tau_{Cs}=$(fit_cs.uncs[5])\$ ms", color="C1")
title("Coprop Rabi flopping")

NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
