#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

Params1 = [0:5:100;]
Params2 = [5:5:100;]
Params3 = [7:7:140;]
Params4 = [10:10:200;]
Params5 = [15:15:300;]
Params6 = [15:15:300;]
Params7 = [15:15:300;]
Params8 = [10:10:200;]
Params9 = [15:15:300;]

offset1 = length(Params1)
offset2 = offset1 + length(Params2)
offset3 = offset2 + length(Params3)
offset4 = offset3 + length(Params4)
offset5 = offset4 + length(Params5)
offset6 = offset5 + length(Params6)
offset7 = offset6 + length(Params7)
offset8 = offset7 + length(Params8)
offset9 = offset8 + length(Params9)

Idx1 = [1:offset1;]
Idx2 = [1; (offset1 + 1):offset2]
Idx3 = [1; (offset2 + 1):offset3]
Idx4 = [1; (offset3 + 1):offset4]
Idx5 = [1; (offset4 + 1):offset5]
Idx6 = [1; (offset5 + 1):offset6]
Idx7 = [1; (offset6 + 1):offset7]
Idx8 = [1; (offset7 + 1):offset8]
Idx9 = [1; (offset8 + 1):offset9]

unshift!(Params2, 0)
unshift!(Params3, 0)
unshift!(Params4, 0)
unshift!(Params5, 0)
unshift!(Params6, 0)
unshift!(Params7, 0)
unshift!(Params8, 0)
unshift!(Params9, 0)

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
plot_params(Params1, Idx1, fmt="o-", label="axial -2")
plot_params(Params2, Idx2, fmt="o-", label="axial -3")
grid()
legend()
ylim([0, ylim()[2]])
title("Axial -2/-3")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_a1_-23.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_-23.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params3, Idx3, fmt="o-", label="axial -4")
plot_params(Params4, Idx4, fmt="o-", label="axial -5")
grid()
legend()
ylim([0, 0.4])
title("Axial -4/-5")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_a1_-45.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_-45.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params5, Idx5, fmt="o-", label="axial -6")
plot_params(Params6, Idx6, fmt="o-", label="axial -7")
plot_params(Params7, Idx7, fmt="o-", label="axial -8")
grid()
legend()
ylim([0, 0.3])
title("Axial -6/-7/-8")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_a1_-678.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_-678.svg", bbox_inches="tight", transparent=true)
close()

figure()
plot_params(Params8, Idx8, fmt="o-", label="axial -1")
plot_params(Params9, Idx9, fmt="o-", label="axial -2")
grid()
legend()
ylim([0, ylim()[2]])
title("Axial -1/-2 (narrow)")
xlabel("t/\$\\mu s\$")
savefig("$(prefix)_a1_-12n.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_-12n.svg", bbox_inches="tight", transparent=true)
close()

# show()
