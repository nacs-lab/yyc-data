#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

Params1 = [linspace(-18.93, -18.71, 12); linspace(-18.125, -17.425, 29)]
Params2 = [linspace(-18.87, -19.07, 11); linspace(-17.925, -17.1, 34)]
Params3 = [linspace(-18.93, -18.71, 12); linspace(-18.125, -17.425, 29)]
Params4 = [linspace(-18.87, -19.07, 11); linspace(-17.925, -17.1, 34)]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)
offset3 = offset2 + length(Params3)
offset4 = offset3 + length(Params4)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]
Idx3 = [(offset2 + 1):offset3;]
Idx4 = [(offset3 + 1):offset4;]

function plot_params(Params, Idx; kws...)
    perm = sortperm(Params)
    Params = Params[perm]
    Idx = Idx[perm]
    Ratios = ratios[Idx, 2]
    Uncs = uncs[Idx, 2]
    errorbar(Params, Ratios, Uncs; kws...)
end

const prefix = ARGS[2]

figure()
plot_params(Params3, Idx3, fmt="ro-", label="Before")
plot_params(Params1, Idx1, fmt="bo-", label="After")
grid()
legend()
ylim([0, 1])
title("Radial 2")
xlabel("\$\\delta\$/MHz")
savefig("$(prefix)_r2.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r2.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params4, Idx4, fmt="ro-", label="Before")
plot_params(Params2, Idx2, fmt="bo-", label="After")
grid()
legend()
ylim([0, 1])
title("Radial 3")
xlabel("\$\\delta\$/MHz")
savefig("$(prefix)_r3.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r3.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params4, Idx4, fmt="ro-")
grid()
ylim([0, 1])
xlabel("\$\\delta\$/MHz")
savefig("$(prefix)_r3_before.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r3_before.svg", bbox_inches="tight", transparent=true)
close()
