#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181116_225125.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const data_nacs_a = NaCsData.select_count(params_a, logicals_a,
                                          NaCsData.select_single((1, 2), (3, 4)))
const times_a = Float64[0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]
const ntimes_a = length(times_a)
const freqs_a = 298.121:0.001:298.139
const spec_a = (1.0:(ntimes_a * length(freqs_a)), 0.0:0.0)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

const iname_b = joinpath(@__DIR__, "data", "data_20181117_102629.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const data_nacs_b = NaCsData.select_count(params_b, logicals_b,
                                          NaCsData.select_single((1, 2), (3, 4)))
const times_b = Float64[0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]
const ntimes_b = length(times_b)
const freqs_b = 298.126:0.0005:298.134
const spec_b = (1.0:(ntimes_b * length(freqs_b)), 0.0:0.0)
const split_nacs_b = NaCsData.split_data(data_nacs_b, spec_b)

const iname_c = joinpath(@__DIR__, "data", "data_20181118_000709.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const data_nacs_c = NaCsData.select_count(params_c, logicals_c,
                                          NaCsData.select_single((1, 2), (3, 4)))
const times_c = Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
const ntimes_c = length(times_c)
const freqs_c = 298.126:0.0005:298.134
const spec_c = (1.0:(ntimes_c * length(freqs_c)), 0.0:0.0)
const split_nacs_c = NaCsData.split_data(data_nacs_c, spec_c)


const time_datas_a =
    [[split_nacs_a[2];
      NaCsData.map_params((i, v)->times_a[i],
                          split_nacs_a[1][((i - 1) * ntimes_a + 1):((i - 1) * ntimes_a + ntimes_a)])]
     for i in 1:length(freqs_a)]
const freq_datas_a =
    [NaCsData.map_params((i, v)->freqs_a[i], split_nacs_a[1][i:ntimes_a:end])
     for i in 1:ntimes_a]

const time_datas_b =
    [[split_nacs_b[2];
      NaCsData.map_params((i, v)->times_b[i],
                          split_nacs_b[1][((i - 1) * ntimes_b + 1):((i - 1) * ntimes_b + ntimes_b)])]
     for i in 1:length(freqs_b)]
const freq_datas_b =
    [NaCsData.map_params((i, v)->freqs_b[i], split_nacs_b[1][i:ntimes_b:end])
     for i in 1:ntimes_b]

const time_datas_c =
    [[split_nacs_c[2];
      NaCsData.map_params((i, v)->times_c[i],
                          split_nacs_c[1][((i - 1) * ntimes_c + 1):((i - 1) * ntimes_c + ntimes_c)])]
     for i in 1:length(freqs_c)]
const freq_datas_c =
    [NaCsData.map_params((i, v)->freqs_c[i], split_nacs_c[1][i:ntimes_c:end])
     for i in 1:ntimes_c]

const time_datas = copy(time_datas_b)
const freq_datas = copy(freq_datas_b)

for i in 1:length(times_b)
    @assert times_b[i] == times_a[i]
    freq_datas[i] = [freq_datas[i]; freq_datas_a[i]]
end

for i in 1:length(freqs_b)
    fb = freqs_b[i]
    j = findfirst(x->x == fb, freqs_a)
    if j === nothing || j === 0
        continue
    end
    time_datas[i] = [time_datas[i]; time_datas_a[j]]
end

const prefix = joinpath(@__DIR__, "imgs", "data_20181116_raman_2d_300")

figure()
for data in time_datas
    NaCsPlot.plot_survival_data(data, fmt=".-")
end
grid()
ylim([0, 1])
title("Raman time")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time")

figure()
for data in freq_datas
    NaCsPlot.plot_survival_data(data, fmt=".-")
end
grid()
ylim([0, 1])
title("Raman spectrum")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_freq")

figure()
for data in time_datas_c
    NaCsPlot.plot_survival_data(data, fmt=".-")
end
grid()
ylim([0, 1])
title("Raman time")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time_short")

figure()
for data in freq_datas_c
    NaCsPlot.plot_survival_data(data, fmt=".-")
end
grid()
ylim([0, 1])
title("Raman spectrum")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_freq_short")

figure()
for i in 1:length(time_datas)
    if freqs_b[i] != 298.130
        continue
    end
    global data_130
    data_130 = time_datas[i]
end
for i in 1:length(time_datas_c)
    if freqs_c[i] != 298.130
        continue
    end
    global data_130
    data_130 = [data_130; time_datas_c[i]]
end
NaCsPlot.plot_survival_data(data_130, fmt=".-")
grid()
ylim([0, 1])
title("Raman time @ 298.13MHz")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time_130khz")

figure()
data_1ms = NaCsData.map_params((i, v)->(v - 298) * 1000, [freq_datas[2]; freq_datas_c[end]])
NaCsPlot.plot_survival_data(data_1ms, fmt=".-")
grid()
ylim([0, 1])
title("Raman spectrum @ 1ms")
xlabel("Detuning (298XXX kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_freq_1ms")

NaCsPlot.maybe_show()
