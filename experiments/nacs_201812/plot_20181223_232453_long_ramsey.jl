#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181223_232453.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [(0:0.5:25,
                   (0:0.5:25) .+ 50,
                   (0:0.5:25) .+ 50 * 2,
                   (0:0.5:25) .+ 50 * 3,
                   (0:0.5:25) .+ 50 * 4,
                   (0:0.5:25) .+ 50 * 6,
                   (0:0.5:50) .+ 50 * 8)]
const specs_cs = [(0:0.016:0.8,
                   (0:0.016:0.8) .+ 2,
                   (0:0.016:0.8) .+ 4,
                   (0:0.016:0.8) .+ 6,
                   (0:0.016:0.8) .+ 9,
                   (0:0.016:0.8) .+ 12,
                   0:8.0:800)]
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

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na = datas_na[1]
data_cs_diff = datas_cs[1][1:6]
data_cs_ax = datas_cs[1][7]

model_sin(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4])
model_sinexp(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4]) .* exp.(.-x ./ p[5])

fit_na = fit_survival(model_sin, data_na[1], [0.4, 0.2, 0.5, -1.7])
fit_na = fit_survival(model_sin, [data_na[1]; data_na[2]], fit_na.param)
fit_na = fit_survival(model_sinexp, [data_na[1]; data_na[2]; data_na[3]], [fit_na.param; 200])
fit_na = fit_survival(model_sinexp, [data_na[1]; data_na[2]; data_na[3]; data_na[4];
                                     data_na[5]; data_na[6]; data_na[7]], fit_na.param)

fit_cs = fit_survival(model_sin, data_cs_diff[1], [0.55, 0.25, 6.3, 1.5])
fit_cs = fit_survival(model_sin, [data_cs_diff[1]; data_cs_diff[2]], fit_cs.param)
fit_cs = fit_survival(model_sinexp, [data_cs_diff[1]; data_cs_diff[2]; data_cs_diff[3]],
                      [fit_cs.param; 7])
fit_cs = fit_survival(model_sinexp, [data_cs_diff[1]; data_cs_diff[2]; data_cs_diff[3];
                                     data_cs_diff[4]; data_cs_diff[5]; data_cs_diff[6]],
                      fit_cs.param)

const prefix = joinpath(@__DIR__, "imgs", "data_20181223_232453_long_ramsey")

# Na

f, axs = subplots(1, 5, sharey=true, facecolor="w", figsize=(16, 4.8))
subplots_adjust(wspace=0.05, hspace=0.05)

d = 0.015 # how big to make the diagonal lines in axes coordinates
sca(axs[1])
axs[1][:spines]["right"][:set_visible](false)
axs[1][:yaxis][:tick_left]()
axs[1][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_na[1], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
grid()
xlim([0, 26])
kwargs = Dict(:transform=>axs[1][:transAxes], :color=>"k", :clip_on=>false)
axs[1][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[1][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)
ylabel("Survival")

sca(axs[2])
axs[2][:spines]["left"][:set_visible](false)
axs[2][:spines]["right"][:set_visible](false)
axs[2][:yaxis][:set_ticks_position]("none")
axs[2][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_na[2], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
grid()
xlim([49.5, 75.5])
kwargs = Dict(:transform=>axs[2][:transAxes], :color=>"k", :clip_on=>false)
axs[2][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[2][:plot]((-d, d), (-d, d); kwargs...)
axs[2][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[2][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[3])
axs[3][:spines]["left"][:set_visible](false)
axs[3][:spines]["right"][:set_visible](false)
axs[3][:yaxis][:set_ticks_position]("none")
axs[3][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_na[3], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
grid()
xlim([99.5, 125.5])
kwargs = Dict(:transform=>axs[3][:transAxes], :color=>"k", :clip_on=>false)
axs[3][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[3][:plot]((-d, d), (-d, d); kwargs...)
axs[3][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[3][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[4])
axs[4][:spines]["left"][:set_visible](false)
axs[4][:spines]["right"][:set_visible](false)
axs[4][:yaxis][:set_ticks_position]("none")
axs[4][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_na[4], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
grid()
xlim([149.5, 175.5])
kwargs = Dict(:transform=>axs[4][:transAxes], :color=>"k", :clip_on=>false)
axs[4][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[4][:plot]((-d, d), (-d, d); kwargs...)
axs[4][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[4][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

# sca(axs[5])
# axs[5][:spines]["left"][:set_visible](false)
# axs[5][:spines]["right"][:set_visible](false)
# axs[5][:yaxis][:set_ticks_position]("none")
# axs[5][:tick_params](labelright="off")
# NaCsPlot.plot_survival_data(data_na[5], fmt="C0.")
# plot(fit_na.plotx, fit_na.ploty, "C0-")
# grid()
# xlim([199.5, 225.5])
# kwargs = Dict(:transform=>axs[5][:transAxes], :color=>"k", :clip_on=>false)
# axs[5][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
# axs[5][:plot]((-d, d), (-d, d); kwargs...)
# axs[5][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
# axs[5][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

# sca(axs[6])
# axs[6][:spines]["left"][:set_visible](false)
# axs[6][:spines]["right"][:set_visible](false)
# axs[6][:yaxis][:set_ticks_position]("none")
# axs[6][:tick_params](labelright="off")
# NaCsPlot.plot_survival_data(data_na[6], fmt="C0.")
# plot(fit_na.plotx, fit_na.ploty, "C0-")
# grid()
# xlim([299.5, 325.5])
# kwargs = Dict(:transform=>axs[6][:transAxes], :color=>"k", :clip_on=>false)
# axs[6][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
# axs[6][:plot]((-d, d), (-d, d); kwargs...)
# axs[6][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
# axs[6][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[5])
axs[5][:spines]["left"][:set_visible](false)
axs[5][:yaxis][:tick_right]()
NaCsPlot.plot_survival_data(data_na[5], fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
grid()
xlim([199.5, 225.5])
kwargs = Dict(:transform=>axs[5][:transAxes], :color=>"k", :clip_on=>false)
axs[5][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[5][:plot]((-d, d), (-d, d); kwargs...)

f[:add_subplot](111, frameon=false)
tick_params(labelcolor="none", top=false, bottom=false, left=false, right=false)
xlabel("Raman time (\$\\mu\$s)")
xlim([0, 1])
ylim([0, 1])
title("Na Y motional-only")

text(0.6, 0.75, "\$\\nu=$(fit_na.uncs[3] * 1000)\$kHz\n" *
     "\$\\tau=$(fit_na.uncs[5])\\mu\$s", color="C0")

NaCsPlot.maybe_save("$(prefix)_na")

# Cs

f, axs = subplots(1, 6, sharey=true, facecolor="w", figsize=(16, 4.8))
subplots_adjust(wspace=0.05, hspace=0.05)

d = 0.015 # how big to make the diagonal lines in axes coordinates
sca(axs[1])
axs[1][:spines]["right"][:set_visible](false)
axs[1][:yaxis][:tick_left]()
axs[1][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_cs_diff[1], fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
xlim([0, 0.9])
kwargs = Dict(:transform=>axs[1][:transAxes], :color=>"k", :clip_on=>false)
axs[1][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[1][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)
ylabel("Survival")

sca(axs[2])
axs[2][:spines]["left"][:set_visible](false)
axs[2][:spines]["right"][:set_visible](false)
axs[2][:yaxis][:set_ticks_position]("none")
axs[2][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_cs_diff[2], fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
xlim([1.95, 2.85])
kwargs = Dict(:transform=>axs[2][:transAxes], :color=>"k", :clip_on=>false)
axs[2][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[2][:plot]((-d, d), (-d, d); kwargs...)
axs[2][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[2][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[3])
axs[3][:spines]["left"][:set_visible](false)
axs[3][:spines]["right"][:set_visible](false)
axs[3][:yaxis][:set_ticks_position]("none")
axs[3][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_cs_diff[3], fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
xlim([3.95, 4.85])
kwargs = Dict(:transform=>axs[3][:transAxes], :color=>"k", :clip_on=>false)
axs[3][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[3][:plot]((-d, d), (-d, d); kwargs...)
axs[3][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[3][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[4])
axs[4][:spines]["left"][:set_visible](false)
axs[4][:spines]["right"][:set_visible](false)
axs[4][:yaxis][:set_ticks_position]("none")
axs[4][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_cs_diff[4], fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
xlim([5.95, 6.85])
kwargs = Dict(:transform=>axs[4][:transAxes], :color=>"k", :clip_on=>false)
axs[4][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[4][:plot]((-d, d), (-d, d); kwargs...)
axs[4][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[4][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[5])
axs[5][:spines]["left"][:set_visible](false)
axs[5][:spines]["right"][:set_visible](false)
axs[5][:yaxis][:set_ticks_position]("none")
axs[5][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_cs_diff[5], fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
xlim([8.95, 9.85])
kwargs = Dict(:transform=>axs[5][:transAxes], :color=>"k", :clip_on=>false)
axs[5][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[5][:plot]((-d, d), (-d, d); kwargs...)
axs[5][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[5][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[6])
axs[6][:spines]["left"][:set_visible](false)
axs[6][:yaxis][:tick_right]()
NaCsPlot.plot_survival_data(data_cs_diff[6], fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
xlim([11.95, 12.85])
kwargs = Dict(:transform=>axs[6][:transAxes], :color=>"k", :clip_on=>false)
axs[6][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[6][:plot]((-d, d), (-d, d); kwargs...)

f[:add_subplot](111, frameon=false)
tick_params(labelcolor="none", top=false, bottom=false, left=false, right=false)
xlabel("Raman time (ms)")
xlim([0, 1])
ylim([0, 1])
title("Cs radial difference")

text(0.6, 0.05, "\$\\nu=$(fit_cs.uncs[3])\$kHz\n" *
     "\$\\tau=$(fit_cs.uncs[5])\$ms", color="C0")

NaCsPlot.maybe_save("$(prefix)_cs_diff")

figure()
NaCsPlot.plot_survival_data(data_cs_ax, fmt="C0")
grid()
xlim([0, 850])
title("Cs axial-only Ramsey")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ax")

NaCsPlot.maybe_show()
