#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200401_182808.mat",
                "data_20200401_220309.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [[0, 0.08, 0.16, 0.24, 0.3] .* 2,
               [0, 0.08, 0.16, 0.24, 0.3] .* 2]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function merge_unc(avg, unc)
    D = sum(unc.^-2)
    N = sum(avg .* unc.^-2)
    return N / D, 1 / sqrt(D)
end

function NaCsCalc.mean(uncs::AbstractArray{T}) where T<:Unc
    return Unc(merge_unc([u.a for u in uncs],
                         [u.s for u in uncs])...)
end

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

const prefix = joinpath(@__DIR__, "imgs", "data_20200401_182808_molecule_lifetime_3322")

function model_expoff(x, p)
    p[3] .+ p[1] .* exp.(.- x .* p[2])
end
function model_exp(x, p)
    p[1] .* exp.(.- x .* p[2])
end
fit1 = fit_survival(model_expoff, datas_nacs[1], [0.1, 6, 0.01])
fit2 = fit_survival(model_expoff, datas_nacs[2], [0.1, 6, 0.01])

# @show fit1.uncs
# @show fit2.uncs

const gamma_avg = mean([fit1.uncs[2], fit2.uncs[2]])

figure()
NaCsPlot.plot_survival_data(datas_nacs[1], fmt="C0s")
plot(fit1.plotx, fit1.ploty, "C0-")
NaCsPlot.plot_survival_data(datas_nacs[2], fmt="C1s")
plot(fit2.plotx, fit2.ploty, "C1-")
xlim([0, 0.63])
ylim([0.035, 0.14])
text(0.25, 0.08, "\$\\gamma=2\\pi\\times$(fit1.uncs[2] / 2π)\$ kHz", color="C0")
text(0.25, 0.07, "\$\\gamma=2\\pi\\times$(fit2.uncs[2] / 2π)\$ kHz", color="C1")
text(0.18, 0.098, "\$\\gamma_{avg}=2\\pi\\times$(gamma_avg / 2π)\$ kHz")
title("306560 GHz, 6 mW")
grid()
xlabel("Molecule Hold Time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
