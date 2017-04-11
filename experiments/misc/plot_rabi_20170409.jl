#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
using DataStructures
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const iname_a = joinpath(@__DIR__, "data", "data_20170404_133229.csv")
const iname_b = joinpath(@__DIR__, "data", "data_20170405_183620.csv")
const iname_c = joinpath(@__DIR__, "data", "data_20170409_141913.csv")
const iname_d = joinpath(@__DIR__, "data", "data_20170409_162137.csv")
const iname_e = joinpath(@__DIR__, "data", "data_20170409_184250.csv")
const iname_f = joinpath(@__DIR__, "data", "data_20170409_233639.csv")

# Before Carrier
data_a = NaCsData.load_count_csv(iname_a)
# Before +1
data_b = NaCsData.load_count_csv(iname_b)
# After a1 +1
data_c = NaCsData.load_count_csv(iname_c)
# After a1 0
data_d = NaCsData.load_count_csv(iname_d)
# After r2 0, +1 (<0 is carrier)
data_e = NaCsData.load_count_csv(iname_e)
# After r3 0, +1
data_f = NaCsData.load_count_csv(iname_f)

# Before Carrier
spec_a = OrderedDict(
    :unused=>(linspace(-17.53, -17.745, 11),
              linspace(-17.14, -17.44, 11),
              linspace(-18.40, -17.90, 51),
              linspace(-18.40, -18.35, 11)),
    :after=>(0:2:100, 2:2:100, 10:10:400),
    :before=>(2:2:30, 2:2:30, 10:10:140)
)
# Before +1
spec_b = OrderedDict(
    :after=>(0:7:280, 8:8:240, 25:25:450),
    :before=>(7:7:140, 8:8:160, 25:25:400)
)
# After r3 0, +1
spec_f = (0:6:180, 2:2:80)

data_a = NaCsData.split_data(data_a, spec_a)
data_b = NaCsData.split_data(data_b, spec_b)
data_f = NaCsData.split_data(data_f, spec_f)

function plot_data(data, scale=1; kws...)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = ratios[perm, 2] .* scale
    uncs = uncs[perm, 2] .* scale
    errorbar(params, ratios, uncs; kws...)
end

const save_fig = get(ENV, "NACS_SAVE_FIG", "true") == "true"

function maybe_save(name)
    if save_fig
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if !save_fig
        show()
    end
end

const prefix = joinpath(@__DIR__, "imgs", "data_rabi_20170409")

data_before_r2_0 = [data_a[:after][1][1]; data_a[:before][1]]
data_before_r3_0 = [data_a[:after][1][1]; data_a[:before][2]]
data_before_a1_0 = [data_a[:after][1][1]; data_a[:before][3]]

data_before_r2_p1 = [data_b[:after][1][1]; data_b[:before][1]]
data_before_r3_p1 = [data_b[:after][1][1]; data_b[:before][2]]
data_before_a1_p1 = [data_b[:after][1][1]; data_b[:before][3]]

data_after_r2_0 = NaCsData.map_params((i, v)->-v * 1e6, data_e[data_e.params .<= 0])
data_after_r2_p1 = NaCsData.map_params((i, v)->v * 1e6, data_e[data_e.params .>= 0])

data_after_r3_0 = [data_f[1][1]; data_f[2]]
data_after_r3_p1 = data_f[1]

data_after_a1_0 = NaCsData.map_params((i, v)->v * 1e6, data_d)
data_after_a1_p1 = NaCsData.map_params((i, v)->v * 1e6, data_c)

figure()
plot_data(data_after_r2_0, 1 / 0.85, fmt="bo-", label="After")
plot_data(data_before_r2_0, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 2 carrier")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r2_0")

figure()
plot_data(data_after_r3_0, 1 / 0.85, fmt="bo-", label="After")
plot_data(data_before_r3_0, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 3 carrier")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r3_0")

figure()
plot_data(data_after_a1_0, 1 / 0.85, fmt="bo-", label="After")
plot_data(data_before_a1_0, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Axial carrier")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_a1_0")

figure()
plot_data(data_after_r2_p1, 1 / 0.85, fmt="bo-", label="After")
plot_data(data_before_r2_p1, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 2 heating")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r2_p1")

figure()
plot_data(data_after_r3_p1, 1 / 0.85, fmt="bo-", label="After")
plot_data(data_before_r3_p1, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 3 heating")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r3_p1")

figure()
plot_data(data_after_a1_p1, 1 / 0.85, fmt="bo-", label="After")
plot_data(data_before_a1_p1, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Axial heating")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_a1_p1")

maybe_show()
