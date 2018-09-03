#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180810_150954.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a,
                                  NaCsData.select_single((1, -2), (3,)))
data_nacs_a = NaCsData.select_count(params_a, logicals_a,
                                    NaCsData.select_single((1, 2), (3,)))

const spec_a = OrderedDict(
    :n0=>-70.0:-11.0,
    :n1=>-70.0:-11.0,
)

const split_na_a = NaCsData.split_data(data_na_a, spec_a)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180810_150954_interaction_shift_n1")

data_na_n0 = split_na_a[:n0]
data_na_n1 = split_na_a[:n1]
data_nacs_n0 = split_nacs_a[:n0]
data_nacs_n1 = split_nacs_a[:n1]

const plt_data_dir = joinpath(@__DIR__, "plot_data")
mkpath(plt_data_dir, 0o755)
const plt_data_prefix = joinpath(plt_data_dir, "data_20180810_150954_interaction_shift_n1")

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

write_datacsv("$(plt_data_prefix)_na_n0", data_na_n0)
write_datacsv("$(plt_data_prefix)_nacs_na_n0", data_nacs_n0)
write_datacsv("$(plt_data_prefix)_na_n1", data_na_n1)
write_datacsv("$(plt_data_prefix)_nacs_na_n1", data_nacs_n1)

figure()
NaCsPlot.plot_survival_data(data_na_n0, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_nacs_n0, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("n=0 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_n0")

figure()
NaCsPlot.plot_survival_data(data_na_n1, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_nacs_n1, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("n=1 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_n1")

figure()
NaCsPlot.plot_survival_data(data_na_n0, fmt="C0.-", label="n=0")
NaCsPlot.plot_survival_data(data_na_n1, fmt="C1.-", label="n=1")
grid()
legend()
ylim([0, 0.8])
title("Na only")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_nacs_n0, fmt="C0.-", label="n=0")
NaCsPlot.plot_survival_data(data_nacs_n1, fmt="C1.-", label="n=1")
grid()
legend()
ylim([0, 0.8])
title("Na + Cs")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs")

NaCsPlot.maybe_show()
