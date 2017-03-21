#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

errorbar(params * 1e-6, ratios[:, 2], uncs[:, 2], fmt="o-")
ylim([0, ylim()[2]])
grid()
xlabel("f/MHz")

axvline(-18.430, color="r")
axvline(-18.364, color="r")
axvline(-18.297, color="r")
axvline(-18.230, color="r")
axvline(-18.163, color="r")
axvline(-18.097, color="r")
axvline(-18.030, color="r")
axvline(-17.963, color="r")

const prefix = ARGS[2]
savefig("$(prefix).png", bbox_inches="tight", transparent=true)
savefig("$(prefix).svg", bbox_inches="tight", transparent=true)
close()
