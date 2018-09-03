#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180808_235802.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180809_093529.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20180809_103153.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const iname_d = joinpath(@__DIR__, "data", "data_20180809_120431.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)
const iname_e = joinpath(@__DIR__, "data", "data_20180809_135451.mat")
const params_e, logicals_e = NaCsData.load_striped_mat(iname_e)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, -2), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((-1, 2), (4,)))
data_na2_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3,)))
data_cs2_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (4,)))

data_na_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, -2), (3,)))
data_cs_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((-1, 2), (4,)))
data_na2_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3,)))
data_cs2_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (4,)))

data_na_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, -2), (3,)))
data_cs_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((-1, 2), (4,)))
data_na2_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (3,)))
data_cs2_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (4,)))

data_na_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, -2), (3,)))
data_cs_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((-1, 2), (4,)))
data_na2_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, 2), (3,)))
data_cs2_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, 2), (4,)))

data_na_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, -2), (3,)))
data_cs_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((-1, 2), (4,)))
data_na2_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, 2), (3,)))
data_cs2_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, 2), (4,)))

data_na_a = NaCsData.map_params((i, f)->i, data_na_a)
data_cs_a = NaCsData.map_params((i, f)->i, data_cs_a)
data_na2_a = NaCsData.map_params((i, f)->i, data_na2_a)
data_cs2_a = NaCsData.map_params((i, f)->i, data_cs2_a)

data_na_b = NaCsData.map_params((i, f)->i, data_na_b)
data_cs_b = NaCsData.map_params((i, f)->i, data_cs_b)
data_na2_b = NaCsData.map_params((i, f)->i, data_na2_b)
data_cs2_b = NaCsData.map_params((i, f)->i, data_cs2_b)

const spec_a = -23 .+ linspace(-60, 60, 121)
const spec_b = -23 .+ linspace(-60, 60, 121)
const spec_c = OrderedDict(
    :cs33=>-62.0:-15.0,
    :cs44=>-23 .+ linspace(-50, 50, 101),
)
const spec_d = OrderedDict(
    :cs33=>-62.0:-15.0,
    :cs44=>-23 .+ linspace(-60, 60, 121),
)

split_na_a = NaCsData.split_data(data_na_a, spec_a)
split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
split_na2_a = NaCsData.split_data(data_na2_a, spec_a)
split_cs2_a = NaCsData.split_data(data_cs2_a, spec_a)

split_na_b = NaCsData.split_data(data_na_b, spec_b)
split_cs_b = NaCsData.split_data(data_cs_b, spec_b)
split_na2_b = NaCsData.split_data(data_na2_b, spec_b)
split_cs2_b = NaCsData.split_data(data_cs2_b, spec_b)

split_na_c = NaCsData.split_data(data_na_c, spec_c)
split_cs_c = NaCsData.split_data(data_cs_c, spec_c)
split_na2_c = NaCsData.split_data(data_na2_c, spec_c)
split_cs2_c = NaCsData.split_data(data_cs2_c, spec_c)

split_na_d = NaCsData.split_data(data_na_d, spec_d)
split_cs_d = NaCsData.split_data(data_cs_d, spec_d)
split_na2_d = NaCsData.split_data(data_na2_d, spec_d)
split_cs2_d = NaCsData.split_data(data_cs2_d, spec_d)

split_na_e = NaCsData.split_data(data_na_e, spec_d)
split_cs_e = NaCsData.split_data(data_cs_e, spec_d)
split_na2_e = NaCsData.split_data(data_na2_e, spec_d)
split_cs2_e = NaCsData.split_data(data_cs2_e, spec_d)

const prefix = joinpath(@__DIR__, "imgs", "data_20180808_interaction_shift2")

data_33_na = [split_na_a; split_na_c[:cs33]; split_na_d[:cs33]; split_na_e[:cs33]]
data_33_na2 = [split_na2_a; split_na2_c[:cs33]; split_na2_d[:cs33]; split_na2_e[:cs33]]
data_33_cs = [split_cs_a; split_cs_c[:cs33]; split_cs_d[:cs33]; split_cs_e[:cs33]]
data_33_cs2 = [split_cs2_a; split_cs2_c[:cs33]; split_cs2_d[:cs33]; split_cs2_e[:cs33]]

data_44_na = [split_na_b; split_na_c[:cs44]; split_na_d[:cs44]; split_na_e[:cs44]]
data_44_na2 = [split_na2_b; split_na2_c[:cs44]; split_na2_d[:cs44]; split_na2_e[:cs44]]
data_44_cs = [split_cs_b; split_cs_c[:cs44]; split_cs_d[:cs44]; split_cs_e[:cs44]]
data_44_cs2 = [split_cs2_b; split_cs2_c[:cs44]; split_cs2_d[:cs44]; split_cs2_e[:cs44]]

const plt_data_dir = joinpath(@__DIR__, "plot_data")
mkpath(plt_data_dir, 0o755)
const plt_data_prefix = joinpath(plt_data_dir, "data_20180808_interaction_shift2")

write_datacsv(fname, x, y, err) = open("$(fname).csv", "w") do io
    write(io, "X,Y,Err\n")
    writedlm(io, [x y err], ',')
end

function write_datacsv(fname, data)
    params, _ratios, _uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = _ratios[perm, 2]
    uncs = _uncs[perm, 2]
    write_datacsv(fname, params, ratios, uncs)
end

write_datacsv("$(plt_data_prefix)_33_na", data_33_na)
write_datacsv("$(plt_data_prefix)_33_nacs_na", data_33_na2)

figure()
NaCsPlot.plot_survival_data(data_33_na, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_33_na2, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.9])
title("(3+2 -> 3+1) interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_33")

function plot_diff(na, cs; kws...)
    params_na, _ratios_na, _uncs_na = NaCsData.get_values(na)
    perm_na = sortperm(params_na)
    params_na = params_na[perm_na]
    ratios_na = _ratios_na[perm_na, 2]
    uncs_na = _uncs_na[perm_na, 2]

    params_cs, _ratios_cs, _uncs_cs = NaCsData.get_values(cs)
    perm_cs = sortperm(params_cs)
    params_cs = params_cs[perm_cs]
    ratios_cs = _ratios_cs[perm_cs, 2]
    uncs_cs = _uncs_cs[perm_cs, 2]

    @assert params_na == params_cs

    ratios = ratios_na .+ 0.95 .- ratios_cs
    uncs = sqrt.(uncs_na.^2 .+ uncs_cs.^2)

    errorbar(params_na, ratios, uncs; kws...)
    return params_na, ratios, uncs
end

figure()
diff_44_na = plot_diff(data_44_na, data_44_cs, fmt="C0.-", label="Na only")
diff_44_na2 = plot_diff(data_44_na2, data_44_cs2, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.9])
title("(4+2 -> 4+1) interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_44")

write_datacsv("$(plt_data_prefix)_44_na", diff_44_na...)
write_datacsv("$(plt_data_prefix)_44_nacs_na", diff_44_na2...)

NaCsPlot.maybe_show()
