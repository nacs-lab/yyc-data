#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20171201_002515.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function gen_selector(na, cs, survival_na)
    if survival_na
        @assert na
    else
        @assert cs
    end
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[1] != na || logicals[2] != cs
            return Int[1, 0, 0]
        end
        return Int[1, 1, survival_na ? logicals[3] : logicals[4]]
    end
end

# Only first 9600 sequences are likely to be good
data_na = NaCsData.select_count(params_a, logicals_a, gen_selector(true, false, true), 9600)
data_cs = NaCsData.select_count(params_a, logicals_a, gen_selector(false, true, false), 9600)
data_na_both = NaCsData.select_count(params_a, logicals_a, gen_selector(true, true, true), 9600)
data_cs_both = NaCsData.select_count(params_a, logicals_a, gen_selector(true, true, false), 9600)

const spec_a = OrderedDict(
    :lo_mf=>[0, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000],
    :hi_mf=>[0, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000],
    :mixed=>0.1 * [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000],
    :lo_f=>[0, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000]
)

const split_na = NaCsData.split_data(data_na, spec_a)
const split_cs = NaCsData.split_data(data_cs, spec_a)
const split_na_both = NaCsData.split_data(data_na_both, spec_a)
const split_cs_both = NaCsData.split_data(data_cs_both, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20171201_002515")

figure()
# low stretched
NaCsPlot.plot_survival_data(split_cs[:lo_mf], fmt="C1o-", label="Cs only")
NaCsPlot.plot_survival_data(split_na[:lo_mf], fmt="C0o-", label="Na only")
NaCsPlot.plot_survival_data(split_cs_both[:lo_mf], fmt="C3o-", label="Cs both")
NaCsPlot.plot_survival_data(split_na_both[:lo_mf], fmt="C2o-", label="Na both")
grid()
ylim([0, 1])
legend()
title("Cs (3, 3) and Na (1, 1)")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_lo_mf")

figure()
# high stretched
NaCsPlot.plot_survival_data(split_cs[:hi_mf], fmt="C1o-", label="Cs only")
NaCsPlot.plot_survival_data(split_na[:hi_mf], fmt="C0o-", label="Na only")
NaCsPlot.plot_survival_data(split_cs_both[:hi_mf], fmt="C3o-", label="Cs both")
NaCsPlot.plot_survival_data(split_na_both[:hi_mf], fmt="C2o-", label="Na both")
grid()
ylim([0, 1])
legend()
title("Cs (4, 4) and Na (2, 2)")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_hi_mf")

figure()
# Mixed
NaCsPlot.plot_survival_data(split_cs[:mixed], fmt="C1o-", label="Cs only")
NaCsPlot.plot_survival_data(split_na[:mixed], fmt="C0o-", label="Na only")
NaCsPlot.plot_survival_data(split_cs_both[:mixed], fmt="C3o-", label="Cs both")
NaCsPlot.plot_survival_data(split_na_both[:mixed], fmt="C2o-", label="Na both")
grid()
ylim([0, 1])
legend()
title("Mixed")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_mixed")

figure()
# Low F
NaCsPlot.plot_survival_data(split_cs[:lo_f], fmt="C1o-", label="Cs only")
NaCsPlot.plot_survival_data(split_na[:lo_f], fmt="C0o-", label="Na only")
NaCsPlot.plot_survival_data(split_cs_both[:lo_f], fmt="C3o-", label="Cs both")
NaCsPlot.plot_survival_data(split_na_both[:lo_f], fmt="C2o-", label="Na both")
grid()
ylim([0, 1])
legend()
title("Low F")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_lo_f")

figure()
# Cs only
NaCsPlot.plot_survival_data(split_cs[:lo_mf], fmt="C0o-", label="Low stretch")
NaCsPlot.plot_survival_data(split_cs[:hi_mf], fmt="C1o-", label="High stretch")
NaCsPlot.plot_survival_data(split_cs[:lo_f], fmt="C2o-", label="Low F")
NaCsPlot.plot_survival_data(split_cs[:mixed], fmt="C3o-", label="Mixed")
grid()
ylim([0, 1])
legend()
title("Cs only")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

figure()
# Na only
NaCsPlot.plot_survival_data(split_na[:lo_mf], fmt="C0o-", label="Low stretch")
NaCsPlot.plot_survival_data(split_na[:hi_mf], fmt="C1o-", label="High stretch")
NaCsPlot.plot_survival_data(split_na[:lo_f], fmt="C2o-", label="Low F")
NaCsPlot.plot_survival_data(split_na[:mixed], fmt="C3o-", label="Mixed")
grid()
ylim([0, 1])
legend()
title("Na only")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
# Cs both
NaCsPlot.plot_survival_data(split_cs_both[:lo_mf], fmt="C0o-", label="Low stretch")
NaCsPlot.plot_survival_data(split_cs_both[:hi_mf], fmt="C1o-", label="High stretch")
NaCsPlot.plot_survival_data(split_cs_both[:lo_f], fmt="C2o-", label="Low F")
NaCsPlot.plot_survival_data(split_cs_both[:mixed], fmt="C3o-", label="Mixed")
grid()
ylim([0, 1])
legend()
title("Cs both")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_both")

figure()
# Na both
NaCsPlot.plot_survival_data(split_na_both[:lo_mf], fmt="C0o-", label="Low stretch")
NaCsPlot.plot_survival_data(split_na_both[:hi_mf], fmt="C1o-", label="High stretch")
NaCsPlot.plot_survival_data(split_na_both[:lo_f], fmt="C2o-", label="Low F")
NaCsPlot.plot_survival_data(split_na_both[:mixed], fmt="C3o-", label="Mixed")
grid()
ylim([0, 1])
legend()
title("Na both")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_both")

NaCsPlot.maybe_show()
