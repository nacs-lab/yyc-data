#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const iname = joinpath(@__DIR__, "data", "data_20170404_133229.csv")
const params, ratios, uncs = NaCsData.calc_survival(iname)

Params1 = [linspace(-17.53, -17.745, 11);]
Params2 = [linspace(-17.14, -17.44, 11);]
Params3 = [linspace(-18.40, -17.90, 51);]
Params4 = [linspace(-18.40, -18.35, 11);]
Params5 = [0:2:100;]
Params6 = [2:2:100;]
Params7 = [10:10:400;]
Params8 = [2:2:30;]
Params9 = [2:2:30;]
Params10 = [10:10:140;]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)
offset3 = offset2 + length(Params3)
offset4 = offset3 + length(Params4)
offset5 = offset4 + length(Params5)
offset6 = offset5 + length(Params6)
offset7 = offset6 + length(Params7)
offset8 = offset7 + length(Params8)
offset9 = offset8 + length(Params9)
offset10 = offset9 + length(Params10)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]
Idx3 = [(offset2 + 1):offset3;]
Idx4 = [(offset3 + 1):offset4;]
Idx5 = [(offset4 + 1):offset5;]
Idx6 = [(offset4 + 1); (offset5 + 1):offset6]
Idx7 = [(offset4 + 1); (offset6 + 1):offset7]
Idx8 = [(offset4 + 1); (offset7 + 1):offset8]
Idx9 = [(offset4 + 1); (offset8 + 1):offset9]
Idx10 = [(offset4 + 1); (offset9 + 1):offset10]

unshift!(Params6, 0)
unshift!(Params7, 0)
unshift!(Params8, 0)
unshift!(Params9, 0)
unshift!(Params10, 0)

function plot_params(Params, Idx, scale; kws...)
    perm = sortperm(Params)
    Params = Params[perm]
    Idx = Idx[perm]
    Ratios = ratios[Idx, 2] .* scale
    Uncs = uncs[Idx, 2] .* scale
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

const prefix = joinpath(@__DIR__, "imgs", "data_20170404_133229")

figure()
plot_params(Params5, Idx5, 1 / 0.85, fmt="bo-", label="After")
plot_params(Params8, Idx8, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 2")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r2")

figure()
plot_params(Params6, Idx6, 1 / 0.85, fmt="bo-", label="After")
plot_params(Params9, Idx9, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 3")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r3")

figure()
plot_params(Params7, Idx7, 1 / 0.85, fmt="bo-", label="After")
plot_params(Params10, Idx10, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Axial 1")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_a1")

maybe_show()
