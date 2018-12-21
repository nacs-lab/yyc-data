#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181219_144521.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [OrderedDict(
    :sb_pi2=>0:4.0:400,
    :sb_pi=>0:4.0:400,
    :ca_pi2=>0:4.0:400,
    :ca_pi=>0:4.0:400,
)]
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

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

data_cs_sb_pi2 = datas_cs[1][:sb_pi2]
data_cs_sb_pi = datas_cs[1][:sb_pi]
data_cs_ca_pi2 = datas_cs[1][:ca_pi2]
data_cs_ca_pi = datas_cs[1][:ca_pi]
data_na_sb_pi2 = datas_na[1][:sb_pi2]
data_na_sb_pi = datas_na[1][:sb_pi]
data_na_ca_pi2 = datas_na[1][:ca_pi2]
data_na_ca_pi = datas_na[1][:ca_pi]

const prefix = joinpath(@__DIR__, "imgs", "data_20181219_144521_ramsey_time_order")

figure()
NaCsPlot.plot_survival_data(data_cs_sb_pi2, fmt="C0", label="\$\\pi / 2\$")
NaCsPlot.plot_survival_data(data_cs_sb_pi, fmt="C1", label="\$\\pi\$")
legend()
grid()
ylim([0, 1])
title("Cs X sideband Ramsey")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_sb")

figure()
NaCsPlot.plot_survival_data(data_cs_ca_pi2, fmt="C0", label="\$\\pi / 2\$")
NaCsPlot.plot_survival_data(data_cs_ca_pi, fmt="C1", label="\$\\pi\$")
legend()
grid()
ylim([0, 1])
title("Cs X carrier Ramsey")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ca")

figure()
NaCsPlot.plot_survival_data(data_na_sb_pi2, fmt="C0", label="\$\\pi / 2\$")
NaCsPlot.plot_survival_data(data_na_sb_pi, fmt="C1", label="\$\\pi\$")
legend()
grid()
ylim([0, 1])
title("Na X sideband Ramsey")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_sb")

figure()
NaCsPlot.plot_survival_data(data_na_ca_pi2, fmt="C0", label="\$\\pi / 2\$")
NaCsPlot.plot_survival_data(data_na_ca_pi, fmt="C1", label="\$\\pi\$")
legend()
grid()
ylim([0, 1])
title("Na X carrier Ramsey")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ca")

NaCsPlot.maybe_show()
