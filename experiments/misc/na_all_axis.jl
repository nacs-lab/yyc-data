#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])
radial2 = params .> 0
radial3 = params .< 0
const prefix = ARGS[2]

figure()
errorbar(params[radial3] ./ 1e6, ratios[radial3, 2], uncs[radial3, 2],
         label="Radial 3", fmt="bo-")
legend()
grid()
ylim([0, ylim()[2]])
xlabel("Detuning (MHz)")
savefig("$(prefix)_radial3.png", bbox_inches="tight", transparent=true)

figure()
errorbar(-params[radial2] ./ 1e6, ratios[radial2, 2], uncs[radial2, 2],
         label="Radial2", fmt="bo-")
legend()
grid()
ylim([0, ylim()[2]])
xlabel("Detuning (MHz)")
savefig("$(prefix)_radial2.png", bbox_inches="tight", transparent=true)
# show()
