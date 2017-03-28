#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

Params1 = [linspace(-18.985, -18.785, 11); linspace(-18.15, -17.95, 11)]
Params2 = [linspace(-18.945, -19.145, 11); linspace(-18.00, -17.75, 11)]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]

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
title("Radial 2")
xlabel("\$\\delta\$/MHz")
maybe_save("$(prefix)_r2")

figure()
plot_params(Params2, Idx2, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Radial 3")
xlabel("\$\\delta\$/MHz")
maybe_save("$(prefix)_r3")

maybe_show()
