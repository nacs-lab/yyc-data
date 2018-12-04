#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181129_190943.mat",
                "data_20181130_161126.mat",
                "data_20181130_210243.mat",
                "data_20181201_013202.mat",
                "data_20181201_074818.mat",
                "data_20181201_121903.mat",
                "data_20181202_122209.mat",
                "data_20181202_150927.mat",
                "data_20181202_190105.mat",
                "data_20181202_214335.mat",
                "data_20181202_234931.mat",
                "data_20181203_023107.mat",
                "data_20181203_050649.mat",
                "data_20181203_074151.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 40000,
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = OrderedDict(
    :f092_1=>OrderedDict(:flip=>[16, 32, 64, 128, 256, 512, 1024, 2048],
                         :life=>[16, 32, 64, 128, 256, 512, 1024, 2048]),
    :f092_2=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                         :p30=>[0, 10, 20, 50, 100, 200, 500]),
    :f396=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f456=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f497=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f514=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f573=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f613=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f652=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f678=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f683=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f686=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f689=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]),
    :f692=>OrderedDict(:p20=>[0, 10, 20, 50, 100, 200, 500],
                       :p30=>[0, 10, 20, 50, 100, 200, 500],
                       :flip=>[0, 10, 20, 50, 100, 200, 500],
                       :life=>[0, 10, 20, 50, 100, 200, 500]))

select_datas(datas, selector, maxcnts, specs) =
    OrderedDict(name=>NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
                for (data, maxcnt, (name, spec)) in zip(datas, maxcnts, specs))

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2,), (3, 4,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)
const datas_na_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4)), maxcnts, specs)
const datas_cs_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20181129_pa_rates")

model_exp1(x, p) = p[1] .* exp.(x ./ -p[2])
model_exp2(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3] .* exp.(x ./ -p[4])
model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

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
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

# Na lifetime
for (name, data_na) in datas_na
    global data_na_lifetime
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_2")
        continue
    end
    if @isdefined data_na_lifetime
        data_na_lifetime = [data_na_lifetime; data_na[:life]]
    else
        data_na_lifetime = data_na[:life]
    end
end
fit_na_lifetime = fit_survival(model_exp1, data_na_lifetime, [1.0, 400],
                               use_unc=false, plot_lo=0)
@show Unc.(fit_na_lifetime.param, fit_na_lifetime.unc)
# Na lifetime in 30mW 1038
for (name, data_na_na) in datas_na_na
    global data_na_1038_30_lifetime
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_1")
        continue
    end
    if @isdefined data_na_1038_30_lifetime
        data_na_1038_30_lifetime = [data_na_1038_30_lifetime; data_na_na[:p30]]
    else
        data_na_1038_30_lifetime = data_na_na[:p30]
    end
end
fit_na_1038_30_lifetime = fit_survival(model_exp1, data_na_1038_30_lifetime, [1.0, 4000],
                                       plot_lo=0)
@show Unc.(fit_na_1038_30_lifetime.param, fit_na_1038_30_lifetime.unc)
# Na lifetime in 20mW 1038
for (name, data_na_na) in datas_na_na
    global data_na_1038_20_lifetime
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_1")
        continue
    end
    if @isdefined data_na_1038_20_lifetime
        data_na_1038_20_lifetime = [data_na_1038_20_lifetime; data_na_na[:p20]]
    else
        data_na_1038_20_lifetime = data_na_na[:p20]
    end
end
fit_na_1038_20_lifetime = fit_survival(model_exp1, data_na_1038_20_lifetime, [1.0, 4000],
                                       plot_lo=0)
@show Unc.(fit_na_1038_20_lifetime.param, fit_na_1038_20_lifetime.unc)
# Cs lifetime
for (name, data_cs) in datas_cs
    global data_cs_lifetime
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_2")
        continue
    end
    if @isdefined data_cs_lifetime
        data_cs_lifetime = [data_cs_lifetime; data_cs[:life]]
    else
        data_cs_lifetime = data_cs[:life]
    end
end
for (name, data_cs_cs) in datas_cs_cs
    global data_cs_lifetime
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_1")
        continue
    end
    data_cs_lifetime = [data_cs_lifetime; data_cs_cs[:p30]]
end
fit_cs_lifetime = fit_survival(model_exp1, data_cs_lifetime, [1.0, 4000], plot_lo=0)
@show Unc.(fit_cs_lifetime.param, fit_cs_lifetime.unc)
# Cs lifetime 20mW
for (name, data_cs_cs) in datas_cs_cs
    global data_cs_20_lifetime
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_1")
        continue
    end
    if @isdefined data_cs_20_lifetime
        data_cs_20_lifetime = [data_cs_20_lifetime; data_cs_cs[:p20]]
    else
        data_cs_20_lifetime = data_cs_cs[:p20]
    end
end
fit_cs_20_lifetime = fit_survival(model_exp1, data_cs_20_lifetime, [1.0, 4000], plot_lo=0)
@show Unc.(fit_cs_20_lifetime.param, fit_cs_20_lifetime.unc)
# Cs scattering
for (name, data_cs) in datas_cs
    global data_cs_scatter
    name = String(name)
    @assert startswith(name, "f")
    if endswith(name, "_2")
        continue
    end
    if @isdefined data_cs_scatter
        data_cs_scatter = [data_cs_scatter; data_cs[:flip]]
    else
        data_cs_scatter = data_cs[:flip]
    end
end
# NaCs lifetime

figure()
NaCsPlot.plot_survival_data(data_na_lifetime, fmt="C0.", label="Na 700")
plot(fit_na_lifetime.plotx, fit_na_lifetime.ploty, "C0-")
NaCsPlot.plot_survival_data(data_na_1038_30_lifetime, fmt="C1.", label="Na 1038 30mW")
plot(fit_na_1038_30_lifetime.plotx, fit_na_1038_30_lifetime.ploty, "C1-")
NaCsPlot.plot_survival_data(data_na_1038_20_lifetime, fmt="C2.", label="Na 1038 20mW")
plot(fit_na_1038_20_lifetime.plotx, fit_na_1038_20_lifetime.ploty, "C2-")
NaCsPlot.plot_survival_data(data_cs_lifetime, fmt="C3.", label="Cs 30mW")
plot(fit_cs_lifetime.plotx, fit_cs_lifetime.ploty, "C3-")
NaCsPlot.plot_survival_data(data_cs_20_lifetime, fmt="C4.", label="Cs 20mW")
plot(fit_cs_20_lifetime.plotx, fit_cs_20_lifetime.ploty, "C4-")
grid()
legend(fontsize="small")
ylim([0, 1])
title("Lifetimes")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_lifetimes")

# figure()
# NaCsPlot.plot_survival_data(data_cs_scatter, fmt="C0.-", label="Cs F3")
# grid()
# legend(fontsize="small")
# ylim([0, 1])
# title("Lifetimes")
# xlabel("Time (ms)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_scatter")

# figure()
# for (name, data) in datas_cs
#     name = String(name)
#     @assert startswith(name, "f")
#     if endswith(name, "_2")
#         continue
#     end
#     NaCsPlot.plot_survival_data(data[:flip], fmt=".-", label=name)
#     # if @isdefined data_na_lifetime
#     #     data_na_lifetime = [data_na_lifetime; data_na[:life]]
#     # else
#     #     data_na_lifetime = data_na[:life]
#     # end
# end
# for (name, data) in datas_cs_cs
#     name = String(name)
#     @assert startswith(name, "f")
#     if endswith(name, "_1")
#         continue
#     end
#     NaCsPlot.plot_survival_data(data[:p20], fmt=".-", label=name)
#     # if @isdefined data_na_lifetime
#     #     data_na_lifetime = [data_na_lifetime; data_na[:life]]
#     # else
#     #     data_na_lifetime = data_na[:life]
#     # end
# end
# grid()
# legend()
# ylim([0, 1])
# # title("Cs (4 + 2) -> (3 + 2)")
# # xlabel("Detuning (kHz)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_cs_shift")

figure()
NaCsPlot.plot_survival_data(datas_nacs[:f092_2][:p20], fmt="C0.-")
NaCsPlot.plot_survival_data(datas_nacs[:f092_2][:p30], fmt="C1.-")
NaCsPlot.plot_survival_data(datas_nacs[:f678][:p20], fmt="C2.-")
NaCsPlot.plot_survival_data(datas_nacs[:f678][:p30], fmt="C3.-")
# NaCsPlot.plot_survival_data(datas_nacs[:f692][:p20], fmt="C0.-")
# NaCsPlot.plot_survival_data(datas_nacs[:f692][:p30], fmt="C1.-")
# NaCsPlot.plot_survival_data(datas_nacs[:f689][:p20], fmt="C2.-")
# NaCsPlot.plot_survival_data(datas_nacs[:f689][:p30], fmt="C3.-")
grid()
# legend()
ylim([0, 1])
# title("Cs (4 + 2) -> (3 + 2)")
# xlabel("Detuning (kHz)")
ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_cs_shift")

NaCsPlot.maybe_show()
