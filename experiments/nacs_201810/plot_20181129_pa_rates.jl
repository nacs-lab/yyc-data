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

function get_param_unc(fit)
    unc = estimate_errors(fit)
    return (param=fit.param, unc=unc,
            uncs=Unc.(fit.param, unc, Sci))
end

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
function process_2body_data(data, model, p0)
    fit = fit_survival(model, data, p0, plot_lo=0)
    return (data=data, fit=fit, uncs=Unc.(fit.param, fit.unc))
end
function plot_2body_data(data, color; fmt="o", label=nothing)
    NaCsPlot.plot_survival_data(data.data, fmt=fmt, color=color, label=label)
    plot(data.fit.plotx, data.fit.ploty, "-", color=color)
end
function rate_model(x, p)
    return p[1] .+ p[2] ./ (x .- 696.26).^2
end

data_2body_20 = Dict{Int,Any}()

data_2body_20[92] = process_2body_data(datas_nacs[:f092_2][:p20], model_exp1, [1.0, 500])
data_2body_20[396] = process_2body_data(datas_nacs[:f396][:p20], model_exp1, [1.0, 500])
data_2body_20[456] = process_2body_data(datas_nacs[:f456][:p20], model_exp1, [1.0, 500])
data_2body_20[497] = process_2body_data(datas_nacs[:f497][:p20], model_exp1, [1.0, 500])
data_2body_20[514] = process_2body_data(datas_nacs[:f514][:p20], model_exp1, [1.0, 500])
data_2body_20[573] = process_2body_data(datas_nacs[:f573][:p20], model_exp1, [1.0, 500])
data_2body_20[613] = process_2body_data(datas_nacs[:f613][:p20], model_exp1, [1.0, 500])
data_2body_20[652] = process_2body_data(datas_nacs[:f652][:p20], model_exp1, [1.0, 500])
data_2body_20[678] = process_2body_data(datas_nacs[:f678][:p20], model_exp1, [1.0, 500])
data_2body_20[683] = process_2body_data(datas_nacs[:f683][:p20], model_exp_off, [1.0, 500, 0.2])
data_2body_20[686] = process_2body_data(datas_nacs[:f686][:p20], model_exp_off, [1.0, 500, 0.2])
data_2body_20[689] = process_2body_data(datas_nacs[:f689][:p20], model_exp_off, [1.0, 500, 0.2])
data_2body_20[692] = process_2body_data(datas_nacs[:f692][:p20], model_exp_off, [1.0, 500, 0.2])

freqs_20 = sort!(collect(keys(data_2body_20)))
rates_20 = Float64[]
rate_uncs_20 = Float64[]
for f in freqs_20
    r = 1000 / data_2body_20[f].uncs[2]
    push!(rates_20, r.a)
    push!(rate_uncs_20, r.s)
end
fit_rate_20 = get_param_unc(curve_fit(rate_model, freqs_20, rates_20,
                                      rate_uncs_20.^-(2/3), [2.0, 100.0]))

data_2body_30 = Dict{Int,Any}()

data_2body_30[92] = process_2body_data(datas_nacs[:f092_2][:p30], model_exp1, [1.0, 500])
data_2body_30[396] = process_2body_data(datas_nacs[:f396][:p30], model_exp1, [1.0, 500])
data_2body_30[456] = process_2body_data(datas_nacs[:f456][:p30], model_exp1, [1.0, 500])
data_2body_30[497] = process_2body_data(datas_nacs[:f497][:p30], model_exp1, [1.0, 500])
data_2body_30[514] = process_2body_data(datas_nacs[:f514][:p30], model_exp1, [1.0, 500])
data_2body_30[573] = process_2body_data(datas_nacs[:f573][:p30], model_exp1, [1.0, 500])
data_2body_30[613] = process_2body_data(datas_nacs[:f613][:p30], model_exp1, [1.0, 500])
data_2body_30[652] = process_2body_data(datas_nacs[:f652][:p30], model_exp1, [1.0, 500])
data_2body_30[678] = process_2body_data(datas_nacs[:f678][:p30], model_exp1, [1.0, 500])
data_2body_30[683] = process_2body_data(datas_nacs[:f683][:p30], model_exp_off, [1.0, 500, 0.2])
data_2body_30[686] = process_2body_data(datas_nacs[:f686][:p30], model_exp_off, [1.0, 500, 0.2])
data_2body_30[689] = process_2body_data(datas_nacs[:f689][:p30], model_exp_off, [1.0, 500, 0.2])
data_2body_30[692] = process_2body_data(datas_nacs[:f692][:p30], model_exp_off, [1.0, 500, 0.2])

freqs_30 = sort!(collect(keys(data_2body_30)))
rates_30 = Float64[]
rate_uncs_30 = Float64[]
for f in freqs_30
    r = 1000 / data_2body_30[f].uncs[2]
    push!(rates_30, r.a)
    push!(rate_uncs_30, r.s)
end
fit_rate_30 = get_param_unc(curve_fit(rate_model, freqs_30, rates_30,
                                      rate_uncs_30.^-(2/3), [2.0, 100.0]))

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

figure()
NaCsPlot.plot_survival_data(data_cs_scatter, fmt="C0.-", label="Cs F3")
grid()
legend(fontsize="small")
ylim([0, 1])
title("Lifetimes")
xlabel("Time (ms)")
ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_scatter")

figure()
plot_2body_data(data_2body_20[92], "C0", label="092")
plot_2body_data(data_2body_20[396], "C1", label="396")
plot_2body_data(data_2body_20[456], "C2", label="456")
plot_2body_data(data_2body_20[497], "C3", label="497")
plot_2body_data(data_2body_20[514], "C4", label="514")
plot_2body_data(data_2body_20[573], "C5", label="573")
plot_2body_data(data_2body_20[613], "C6", label="613")
plot_2body_data(data_2body_20[652], "C7", label="652")
plot_2body_data(data_2body_20[678], "C8", label="678")
plot_2body_data(data_2body_20[683], "C9", label="683")
plot_2body_data(data_2body_20[686], "C0", fmt="s", label="686")
plot_2body_data(data_2body_20[689], "C1", fmt="s", label="689")
plot_2body_data(data_2body_20[692], "C2", fmt="s", label="692")
grid()
legend(ncol=4, fontsize="small", labelspacing=0.2, borderpad=0.2,
       handletextpad=0.2, columnspacing=0.3, borderaxespad=0.2)
ylim([0, 1])
title("2 body survival (20mW)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_pa20")

figure()
plot_2body_data(data_2body_30[92], "C0", label="092")
plot_2body_data(data_2body_30[396], "C1", label="396")
plot_2body_data(data_2body_30[456], "C2", label="456")
plot_2body_data(data_2body_30[497], "C3", label="497")
plot_2body_data(data_2body_30[514], "C4", label="514")
plot_2body_data(data_2body_30[573], "C5", label="573")
plot_2body_data(data_2body_30[613], "C6", label="613")
plot_2body_data(data_2body_30[652], "C7", label="652")
plot_2body_data(data_2body_30[678], "C8", label="678")
plot_2body_data(data_2body_30[683], "C9", label="683")
plot_2body_data(data_2body_30[686], "C0", fmt="s", label="686")
plot_2body_data(data_2body_30[689], "C1", fmt="s", label="689")
plot_2body_data(data_2body_30[692], "C2", fmt="s", label="692")
grid()
legend(ncol=4, fontsize="small", labelspacing=0.2, borderpad=0.2,
       handletextpad=0.2, columnspacing=0.3, borderaxespad=0.2)
ylim([0, 1])
title("2 body survival (30mW)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_pa30")

figure()
errorbar(freqs_30, rates_30, rate_uncs_30, fmt="C0o", label="30mW")
plot(linspace(600, 695, 1000), rate_model(linspace(600, 695, 1000), fit_rate_30.param), "C0")
errorbar(freqs_20, rates_20, rate_uncs_20, fmt="C1o", label="20mW")
plot(linspace(600, 695, 1000), rate_model(linspace(600, 695, 1000), fit_rate_20.param), "C1")
text(603, 22, "\$\\Gamma=\\Gamma_0+\\dfrac{a}{(f-696.26)^2}\$\n")
text(603, 5, "\$\\Gamma_0=$(fit_rate_30.uncs[1])\$Hz\n" *
     "\$a=$(fit_rate_30.uncs[2])\$Hz\$\\cdot\$GHz\${}^2\$", color="C0", fontsize="small")
text(603, 14, "\$\\Gamma_0=$(fit_rate_20.uncs[1])\$Hz\n" *
     "\$a=$(fit_rate_20.uncs[2])\$Hz\$\\cdot\$GHz\${}^2\$", color="C1", fontsize="small")
legend(fontsize="small")
xlabel("Tweezer frequency (288XXX GHz)")
ylabel("Decay rate (Hz)")
grid()
xlim([600, 695])
ylim([0, 40])
NaCsPlot.maybe_save("$(prefix)_rates")

NaCsPlot.maybe_show()
