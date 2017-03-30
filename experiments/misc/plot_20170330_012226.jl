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
Params3 = [linspace(-18.55, -18.58, 16); linspace(-18.42, -18.45, 16)]
Params4 = [linspace(10, 50, 11);]
Params5 = [linspace(10, 50, 11);]
Params6 = [linspace(40, 160, 13);]
Params7 = [linspace(-18.55, -18.45, 11); linspace(-12.3, -12.2, 11)]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)
offset3 = offset2 + length(Params3)
offset4 = offset3 + length(Params4)
offset5 = offset4 + length(Params5)
offset6 = offset5 + length(Params6)
offset7 = offset6 + length(Params7)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]
Idx3 = [(offset2 + 1):offset3;]
Idx4 = [(offset3 + 1):offset4;]
Idx5 = [(offset4 + 1):offset5;]
Idx6 = [(offset5 + 1):offset6;]
Idx7 = [(offset6 + 1):offset7;]

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

figure()
plot_params(Params3, Idx3, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Axial 1")
xlabel("\$\\delta\$/MHz")
maybe_save("$(prefix)_a1")

figure()
plot_params(Params4, Idx4, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Radial 2 heating sideband")
xlabel("\$t/\\mu s\$")
maybe_save("$(prefix)_r2_+1t")

figure()
plot_params(Params5, Idx5, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Radial 3 heating sideband")
xlabel("\$t/\\mu s\$")
maybe_save("$(prefix)_r3_+1t")

figure()
plot_params(Params6, Idx6, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Axial 1 heating sideband")
xlabel("\$t/\\mu s\$")
maybe_save("$(prefix)_a1_+1t")

figure()
plot_params(Params7, Idx7, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Co-prop")
xlabel("\$t/\\mu s\$")
maybe_save("$(prefix)_co")

maybe_show()
