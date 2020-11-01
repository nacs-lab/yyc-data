#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using DelimitedFiles

const inames = ["data_20180815_000658.mat",
                "data_20180815_091610.mat",
                "data_20180815_231138.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [
    OrderedDict(
        :n0=>7 .+ linspace(-40, 40, 81),
        :n1=>7 .+ linspace(-40, 40, 81),
    ),
    OrderedDict(
        :n0=>7 .+ linspace(-60.0, 10.0, 81),
        :n1=>7 .+ linspace(-60.0, 10.0, 81),
    ),
    OrderedDict(
        :cs=>linspace(-80.0, 80.0, 161),
        :na=>linspace(-80.0, 80.0, 161),
    ),
]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs_cs = select_datas(datas, NaCsData.select_single((1, 2), (4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (4,)), maxcnts, specs)
const datas_nacs_na = select_datas(datas, NaCsData.select_single((1, 2), (3,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20180815_interaction_shift4")

data_cs_n0 = [datas_cs[1][:n0]; datas_cs[2][:n0]; datas_cs[3][:cs]]
data_cs_n1 = [datas_cs[1][:n1]; datas_cs[2][:n1]]
data_nacs_cs_n0 = [datas_nacs_cs[1][:n0]; datas_nacs_cs[2][:n0]; datas_nacs_cs[3][:cs]]
data_nacs_cs_n1 = [datas_nacs_cs[1][:n1]; datas_nacs_cs[2][:n1]]

data_na_n0 = datas_na[3][:na]
data_nacs_na_n0 = datas_nacs_na[3][:na]

const plt_data_dir = joinpath(@__DIR__, "plot_data")
mkpath(plt_data_dir, mode=0o755)
const plt_data_prefix = joinpath(plt_data_dir, "data_20180815_interaction_shift4")

write_datacsv(fname, data) = open("$(fname).csv", "w") do io
    params, _ratios, _uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = _ratios[perm, 2]
    uncs = _uncs[perm, 2]
    write(io, "X,Y,Err\n")
    writedlm(io, [params ratios uncs], ',')
end

write_datacsv("$(plt_data_prefix)_na_n0", data_na_n0)
write_datacsv("$(plt_data_prefix)_nacs_na_n0", data_nacs_na_n0)
write_datacsv("$(plt_data_prefix)_cs_n0", data_cs_n0)
write_datacsv("$(plt_data_prefix)_nacs_cs_n0", data_nacs_cs_n0)
write_datacsv("$(plt_data_prefix)_cs_n1", data_cs_n1)
write_datacsv("$(plt_data_prefix)_nacs_cs_n1", data_nacs_cs_n1)

figure()
NaCsPlot.plot_survival_data(data_na_n0, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_nacs_na_n0, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("(3+2 -> 3+1) n=0 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_n0")

figure()
NaCsPlot.plot_survival_data(data_cs_n0, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_cs_n0, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) n=0 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_n0")

figure()
NaCsPlot.plot_survival_data(data_cs_n1, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_cs_n1, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) n=1 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_n1")

figure()
NaCsPlot.plot_survival_data(data_cs_n0, fmt="C0.-", label="n=0")
NaCsPlot.plot_survival_data(data_cs_n1, fmt="C1.-", label="n=1")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) Cs only")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

figure()
NaCsPlot.plot_survival_data(data_nacs_cs_n0, fmt="C0.-", label="n=0")
NaCsPlot.plot_survival_data(data_nacs_cs_n1, fmt="C1.-", label="n=1")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) Na + Cs")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_cs")

NaCsPlot.maybe_show()
