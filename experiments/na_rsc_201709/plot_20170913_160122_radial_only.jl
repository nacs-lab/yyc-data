#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using PyPlot
using DataStructures
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const iname_a = joinpath(@__DIR__, "data", "data_20170913_160122.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :radial_y=>(linspace(18, 18.4, 20), linspace(19, 19.4, 20), linspace(19.5, 19.8, 15)),
    :radial_x=>(linspace(18, 18.4, 20), linspace(19, 19.4, 20), linspace(19.5, 19.8, 15)),
)

const split_a = NaCsData.split_data(data_a, spec_a)

function plot_data(data, scale=1; kws...)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = ratios[perm, 2] .* scale
    uncs = uncs[perm, 2] .* scale
    errorbar(params, ratios, uncs; kws...)
end

function maybe_save(name)
    if !interactive()
        dir = dirname(name)
        if !isempty(dir)
            mkpath(dir, 0o755)
        end
        savefig("$name.pdf"; bbox_inches="tight", transparent=true)
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

const prefix = joinpath(@__DIR__, "imgs", "data_20170913_160122")

to_sideband(f) = (i, v)->(v - f) * 1000

data_rx = NaCsData.map_params(to_sideband(18.705), split_a[:radial_x])
data_ry = NaCsData.map_params(to_sideband(18.713), split_a[:radial_y])

figure()
# Without cooling
plot_data(data_rx[1], fmt="C0o-")
plot_data(data_rx[2], fmt="C0o-")
plot_data(data_rx[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
maybe_save("$(prefix)_rx")

figure()
# Without cooling
plot_data(data_ry[1], fmt="C0o-")
plot_data(data_ry[2], fmt="C0o-")
plot_data(data_ry[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
maybe_save("$(prefix)_ry")

maybe_show()
