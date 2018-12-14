#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181213_235613.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [(0:0.4:8,
                   (0:0.4:8) .+ 50,
                   (0:0.4:8) .+ 50 * 2,
                   (0:0.4:8) .+ 50 * 4,
                   (0:0.4:8) .+ 50 * 6,
                   0:0.4:8,
                   (0:0.4:8) .+ 50,
                   (0:0.4:8) .+ 50 * 2,
                   (0:0.4:8) .+ 50 * 4,
                   (0:0.4:8) .+ 50 * 6)]
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

data_na = datas_na[1]
data_na = [[data_na[1]; data_na[6]],
           [data_na[2]; data_na[7]],
           [data_na[3]; data_na[8]],
           [data_na[4]; data_na[9]],
           [data_na[5]; data_na[10]]]

model_sin(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4])
model_sinexp(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4]) .* exp.(.-x ./ p[5])

fit = fit_survival(model_sin, data_na[1], [0.4, 0.3, 0.5, 0.5])
fit = fit_survival(model_sin, [data_na[1]; data_na[2]], fit.param)
fit = fit_survival(model_sinexp, [data_na[1]; data_na[2]; data_na[3]], [fit.param; 200])
fit = fit_survival(model_sinexp, [data_na[1]; data_na[2]; data_na[3];
                                  data_na[4]; data_na[5]], fit.param)

const prefix = joinpath(@__DIR__, "imgs", "data_20181213_235613_na_ramsey")

f, axs = subplots(1, 5, sharey=true, facecolor="w", figsize=(10, 4.8))
subplots_adjust(wspace=0.05, hspace=0.05)

d = 0.015 # how big to make the diagonal lines in axes coordinates
sca(axs[1])
axs[1][:spines]["right"][:set_visible](false)
axs[1][:yaxis][:tick_left]()
axs[1][:tick_params](labelright="off")
NaCsPlot.plot_survival_data(data_na[1], fmt="C0.")
plot(fit.plotx, fit.ploty, "C0-")
grid()
xlim([0, 9])
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
plot(fit.plotx, fit.ploty, "C0-")
grid()
xlim([49, 59])
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
plot(fit.plotx, fit.ploty, "C0-")
grid()
xlim([99, 109])
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
plot(fit.plotx, fit.ploty, "C0-")
grid()
xlim([199, 209])
kwargs = Dict(:transform=>axs[4][:transAxes], :color=>"k", :clip_on=>false)
axs[4][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[4][:plot]((-d, d), (-d, d); kwargs...)
axs[4][:plot]((1 - d, 1 + d), (-d, d); kwargs...)
axs[4][:plot]((1 - d, 1 + d), (1 - d, 1 + d); kwargs...)

sca(axs[5])
axs[5][:spines]["left"][:set_visible](false)
axs[5][:yaxis][:tick_right]()
NaCsPlot.plot_survival_data(data_na[5], fmt="C0.")
plot(fit.plotx, fit.ploty, "C0-")
grid()
xlim([299, 309])
kwargs = Dict(:transform=>axs[5][:transAxes], :color=>"k", :clip_on=>false)
axs[5][:plot]((-d, d), (1 - d, 1 + d); kwargs...)
axs[5][:plot]((-d, d), (-d, d); kwargs...)

f[:add_subplot](111, frameon=false)
tick_params(labelcolor="none", top=false, bottom=false, left=false, right=false)
xlabel("Raman time (us)")
xlim([0, 1])
ylim([0, 1])

text(0.6, 0.7, "\$\\nu=$(fit.uncs[3] * 1000)\$kHz\n" *
     "\$\\tau=$(fit.uncs[5])\\mu\$s", color="C0")

NaCsPlot.maybe_save("$(prefix)")


# text(2, 0.05, "\$\\nu_{579}=$(fit_2.uncs[3] * 1000)\$kHz", color="C1")

NaCsPlot.maybe_show()
