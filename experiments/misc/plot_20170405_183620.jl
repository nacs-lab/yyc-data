#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const iname = joinpath(@__DIR__, "data", "data_20170405_183620.csv")
const params, ratios, uncs = NaCsData.calc_survival(iname)

Params1 = [0:7:280;]
Params2 = [8:8:240;]
Params3 = [25:25:450;]
Params4 = [7:7:140;]
Params5 = [8:8:160;]
Params6 = [25:25:400;]
Params7_1 = [4:4:80;]
Params7_2 = [520:4:600;]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)
offset3 = offset2 + length(Params3)
offset4 = offset3 + length(Params4)
offset5 = offset4 + length(Params5)
offset6 = offset5 + length(Params6)
offset7_1 = offset6 + length(Params7_1)
offset7_2 = offset7_1 + length(Params7_2)

Idx1 = [1:offset1;]
Idx2 = [1; (offset1 + 1):offset2]
Idx3 = [1; (offset2 + 1):offset3]
Idx4 = [1; (offset3 + 1):offset4]
Idx5 = [1; (offset4 + 1):offset5]
Idx6 = [1; (offset5 + 1):offset6]
Idx7_1 = [1; (offset6 + 1):offset7_1]
Idx7_2 = [(offset7_1 + 1):offset7_2;]

unshift!(Params2, 0)
unshift!(Params3, 0)
unshift!(Params4, 0)
unshift!(Params5, 0)
unshift!(Params6, 0)
unshift!(Params7_1, 0)

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

const prefix = joinpath(@__DIR__, "imgs", "data_20170405_183620")

figure()
plot_params(Params1, Idx1, 1 / 0.85, fmt="bo-", label="After")
plot_params(Params4, Idx4, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 2")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r2")

figure()
plot_params(Params2, Idx2, 1 / 0.85, fmt="bo-", label="After")
plot_params(Params5, Idx5, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Radial 3")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r3")

figure()
plot_params(Params3, Idx3, 1 / 0.85, fmt="bo-", label="After")
plot_params(Params6, Idx6, 1 / 0.95, fmt="ro-", label="Before")
grid()
ylim([0, 1.0])
title("Axial 1")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_a1")

figure()
plot_params(Params7_1, Idx7_1, 1 / 0.85, fmt="bo-")
plot_params(Params7_2, Idx7_2, 1 / 0.85, fmt="bo-")
grid()
ylim([0, 1.0])
title("Co-prop")
xlabel("\$t/\\mu s\$")
ylabel("Normalized survival")
maybe_save("$(prefix)_co")

maybe_show()
