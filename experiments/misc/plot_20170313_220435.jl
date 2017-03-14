#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

Params1 = [0:40:400; 1600:40:2000]
offset1 = length(Params1)
Params2 = [linspace(-18.03, -18.08, 11); linspace(-18.942, -18.892, 11)]
offset2 = offset1 + length(Params2)
Params3 = [20:20:200; 800:20:1000]
offset3 = offset2 + length(Params3)
Params4 = [linspace(-17.855, -17.905, 11); linspace(-19.073, -19.123, 11)]
offset4 = offset3 + length(Params4)
Params5 = [20:20:200; 800:20:1000]
offset5 = offset4 + length(Params5)
Params6 = [linspace(-18.54, -18.57, 16); linspace(-18.40, -18.43, 16)]
offset6 = offset5 + length(Params6)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]
Idx3 = [1; (offset2 + 1):offset3]
Idx4 = [(offset3 + 1):offset4;]
Idx5 = [1; (offset4 + 1):offset5]
Idx6 = [(offset5 + 1):offset6;]
unshift!(Params3, 0)
unshift!(Params5, 0)

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
plot_params(Params1, Idx1, fmt="o-")
grid()
ylim([0, 0.8])
title("Radial 2 carrier")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_r2_carrier.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r2_carrier.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params2, Idx2, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("Radial 2 cooling+heating")
xlabel("f/\$MHz\$")
savefig("$(prefix)_r2_+-1.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r2_+-1.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params3, Idx3, fmt="o-")
grid()
ylim([0, 0.8])
title("Radial 3 carrier")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_r3_carrier.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r3_carrier.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params4, Idx4, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("Radial 3 cooling+heating")
xlabel("f/\$MHz\$")
savefig("$(prefix)_r3_+-1.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_r3_+-1.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params5, Idx5, fmt="o-")
grid()
ylim([0, 0.8])
title("Axial 1 carrier")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_a1_carrier.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_carrier.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params6, Idx6, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("Axial 1 cooling+heating")
xlabel("f/\$MHz\$")
savefig("$(prefix)_a1_+-1.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_+-1.svg", bbox_inches="tight", transparent=true)
close()
