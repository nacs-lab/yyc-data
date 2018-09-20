#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

load_mat(name) = NaCsData.load_striped_mat(joinpath(@__DIR__, "data", name))

names = [
    "data_20180919_230048.mat",
    "data_20180920_000758.mat",
    "data_20180920_011532.mat",
    "data_20180920_022920.mat",
    "data_20180920_033832.mat",
    "data_20180920_044858.mat",
    "data_20180920_055809.mat"
]
files = [load_mat(name) for name in names]

split_data(spec, f, ls...) =
    NaCsData.split_data(NaCsData.select_count(f..., NaCsData.select_single(ls...)), spec)

const spec_a = OrderedDict(
    :cool=>[0, 10, 20, 50, 100, 200, 400],
    :nocool=>[0, 10, 20, 50, 100, 200, 400],
)

freqs = 690:696

datas_na = [split_data(spec_a, f, (1, -2), (3, -4)) for f in files]
datas_cs = [split_data(spec_a, f, (-1, 2), (-3, 4)) for f in files]
datas_nacs_na = [split_data(spec_a, f, (1, 2), (3,)) for f in files]
datas_nacs_cs = [split_data(spec_a, f, (1, 2), (4,)) for f in files]

const prefix = joinpath(@__DIR__, "imgs", "data_20180919_2body_lifetimes")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_na[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_cs[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Cs, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_na[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_nocool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_cs[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Cs, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_nocool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_na[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na + Cs, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_na_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_cs[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Na + Cs, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_cs_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_na[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na + Cs, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_na_nocool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_cs[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Na + Cs, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_cs_nocool")

NaCsPlot.maybe_show()
