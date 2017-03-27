#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

Params1 = [linspace(-18.93, -18.71, 12); linspace(-18.125, -17.425, 29)] .- 0.075
Params2 = [linspace(-18.87, -19.07, 11); linspace(-17.925, -17.1, 34)] .- 0.075
Params3 = [linspace(-18.55, -18.58, 16); linspace(-18.42, -18.45, 16)]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)
offset3 = offset2 + length(Params3)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]
Idx3 = [(offset2 + 1):offset3;]

function plot_params(Params, Idx; kws...)
    perm = sortperm(Params)
    Params = Params[perm]
    Idx = Idx[perm]
    Ratios = ratios[Idx, 2]
    Uncs = uncs[Idx, 2]
    errorbar(Params, Ratios, Uncs; kws...)
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

const prefix = ARGS[2]

figure()
plot_params(Params1, Idx1, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Axis 2 (Radial)")
xlabel("\$\\delta\$/MHz")
maybe_save("$(prefix)_r2")

figure()
plot_params(Params2, Idx2, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Axis 3 (Radial)")
xlabel("\$\\delta\$/MHz")
maybe_save("$(prefix)_r3")

figure()
plot_params(Params3, Idx3, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Axis 1 (Axial)")
xlabel("\$\\delta\$/MHz")
maybe_save("$(prefix)_a1")

maybe_show()
