#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190803_125057.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 ]
const specs = [([0, 1, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 5, 10, 20, 50, 100],
                [0, 1, 2, 5, 10, 20, 50, 100])]
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

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

data_nacs = datas_nacs[1]

fits_nacs = [fit_survival(model_exp_off, data, [0.25, 20, 0.05]) for data in data_nacs]

const prefix = joinpath(@__DIR__, "imgs", "data_20190803_125057_lifetime_33_22")

tps = [1, 2, 6, 9, 15];
τs = Float64[]
τs_s = Float64[]

figure()
for i in 1:5
    NaCsPlot.plot_survival_data(data_nacs[i], fmt="C$(i - 1).", label="$(tps[i]) mW")
    plot(fits_nacs[i].plotx, fits_nacs[i].ploty, "C$(i - 1)-")
    push!(τs, fits_nacs[i].param[2])
    push!(τs_s, fits_nacs[i].unc[2])
end
legend(ncol=2, fontsize="small")
grid()
xlim([0, 110])
title("3, 3 + 2, 2 Lifetime")
xlabel("Hold time (ms)")
ylabel("Two body survival")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(tps, τs, τs_s, fmt="C0o-")
grid()
xlim([0, 16])
title("Power dependency of lifetime")
xlabel("Tweezer Power (mW)")
ylabel("3, 3 + 2, 2 Lifetime")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_taus")

NaCsPlot.maybe_show()
