#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit


const inames = ["data_20181205_212646.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [OrderedDict(:r50=>311.3 .+ [-600:100:-400; -330; -270; -220; -180;
                                           -150:20:150;
                                           180; 220; 270; 330; 400:100:600] .* 1e-3,
                           :r200=>311.3 .+ [-600:100:-400; -330; -270; -220; -180;
                                            -150:20:150;
                                            180; 220; 270; 330; 400:100:600] .* 1e-3)]

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

function gen_model(v0, v1)
    dv = v1 - v0
    return (x, p)->v0 .+ dv .* exp.(-p[1] ./ (1 .+ (x .- p[2]).^2 ./ p[3]^2))
end

function model′(x, p)
    return p[4] .- p[1] ./ (1 .+ (x .- p[2]).^2 ./ p[3]^2)
end

fit_r50 = fit_survival(gen_model(0.1, 0.7), datas_nacs[1][:r50], [0.7, 311.3, 0.1])
@show fit_r50.uncs
fit_r200 = fit_survival(gen_model(0.1, 0.7), datas_nacs[1][:r200], [0.7, 311.3, 0.1])
@show fit_r200.uncs

fit_r50′ = fit_survival(model′, datas_nacs[1][:r50], [0.2, 311.3, 0.1, 0.7])
@show fit_r50′.uncs
fit_r200′ = fit_survival(model′, datas_nacs[1][:r200], [0.2, 311.3, 0.1, 0.7])
@show fit_r200′.uncs

# a + b * exp(-r0 / (1 + (f - f0)^2 / Γ^2))

const prefix = joinpath(@__DIR__, "imgs", "data_20181205_212646_raman_ratios")

figure()
NaCsPlot.plot_survival_data(datas_nacs[1][:r50], fmt="C0.", label="0.02")
plot(fit_r50.plotx, fit_r50.ploty, "C0-")
plot(fit_r50′.plotx, fit_r50′.ploty, "C0--")
NaCsPlot.plot_survival_data(datas_nacs[1][:r200], fmt="C1.", label="0.005")
plot(fit_r200.plotx, fit_r200.ploty, "C1-")
plot(fit_r200′.plotx, fit_r200′.ploty, "C1--")
grid()
legend()
xlabel("Raman frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
