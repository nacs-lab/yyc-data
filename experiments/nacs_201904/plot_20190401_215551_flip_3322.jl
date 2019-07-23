#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190401_215551.mat",
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 ]
const specs = [([0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
                [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
                [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
                [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
                [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
                [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5)]
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

const datas_na_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4)), maxcnts, specs)
const datas_cs_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4)), maxcnts, specs)
const datas_nacs_na = select_datas(datas, NaCsData.select_single((1, 2), (3,)), maxcnts, specs)
const datas_nacs_cs = select_datas(datas, NaCsData.select_single((1, 2), (4,)), maxcnts, specs)

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

# 1. 0-shift survival
# 2. shift survival
# 3. 0-shift F=1/F=3
# 4. shift F=1/F=3
# 5. 0-shift (2, 2)/non-(3, 3)
# 6. shift (2, 2)/non-(3, 3)

data_na_lifetime = [datas_na_na[1][1]; datas_na_na[1][2]]
data_na_f1 = [datas_na_na[1][3]; datas_na_na[1][4]]
data_na_f22 = [datas_na_na[1][5]; datas_na_na[1][6]]

data_cs_lifetime0 = datas_cs_cs[1][1]
data_cs_lifetime1 = datas_cs_cs[1][2]
data_cs_f3_0 = datas_cs_cs[1][3]
data_cs_f3_1 = datas_cs_cs[1][4]
data_cs_f44_0 = datas_cs_cs[1][5]
data_cs_f44_1 = datas_cs_cs[1][6]

# data_na2_lifetime0 = datas_nacs_na[1][1]
# data_na2_state0 = datas_nacs_na[1][2]
# data_na2_lifetime1 = datas_nacs_na[1][3]
# data_na2_state1 = datas_nacs_na[1][4]

# data_cs2_lifetime0 = datas_nacs_cs[1][1]
# data_cs2_state0 = datas_nacs_cs[1][2]
# data_cs2_lifetime1 = datas_nacs_cs[1][3]
# data_cs2_state1 = datas_nacs_cs[1][4]

const prefix = joinpath(@__DIR__, "imgs", "data_20190401_215551_flip_3322")

figure()
NaCsPlot.plot_survival_data(data_na_lifetime, fmt="C0.-", label="Na 2,2")
NaCsPlot.plot_survival_data(data_cs_lifetime0, fmt="C1.-", label="Cs 3,3")
NaCsPlot.plot_survival_data(data_cs_lifetime1, fmt="C2.-", label="Cs 4,4")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 5.5])
title("Single-body lifetime")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_single_lifetime")

figure()
NaCsPlot.plot_survival_data(data_na_f1, fmt="C0.-", label="Na 1")
NaCsPlot.plot_survival_data(data_na_f22, fmt="C1.-", label="Na 2,2")
NaCsPlot.plot_survival_data(data_cs_f3_0, fmt="C2.-", label="Cs 3")
NaCsPlot.plot_survival_data(data_cs_f44_0, fmt="C3.-", label="Cs 4,4")
NaCsPlot.plot_survival_data(data_cs_f3_1, fmt="C4.-", label="Cs 3,3")
NaCsPlot.plot_survival_data(data_cs_f44_1, fmt="C5.-", label="Cs 4,4")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 5.5])
title("Single-body spin flip")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_single_state")

# figure()
# NaCsPlot.plot_survival_data(data_na2_lifetime0, fmt="C0.-", label="Na 0")
# NaCsPlot.plot_survival_data(data_na2_lifetime1, fmt="C1.-", label="Na 1")
# NaCsPlot.plot_survival_data(data_cs2_lifetime0, fmt="C2.-", label="Cs 0")
# NaCsPlot.plot_survival_data(data_cs2_lifetime1, fmt="C3.-", label="Cs 1")
# legend(fontsize="small")
# grid()
# ylim([0, 1])
# xlim([0, 5.5])
# title("Two-body lifetime")
# xlabel("Wait time (ms)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_double_lifetime")

# figure()
# NaCsPlot.plot_survival_data(data_na2_state0, fmt="C0.-", label="Na 0")
# NaCsPlot.plot_survival_data(data_na2_state1, fmt="C1.-", label="Na 1")
# NaCsPlot.plot_survival_data(data_cs2_state0, fmt="C2.-", label="Cs 0")
# NaCsPlot.plot_survival_data(data_cs2_state1, fmt="C3.-", label="Cs 1")
# legend(fontsize="small")
# grid()
# ylim([0, 1])
# xlim([0, 5.5])
# title("Two-body spin flip")
# xlabel("Wait time (ms)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_double_state")

NaCsPlot.maybe_show()
