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
xlabel("\$\\delta\$/MHz")

const prefix = ARGS[2]
savefig("$(prefix).png", bbox_inches="tight", transparent=true)
savefig("$(prefix).svg", bbox_inches="tight", transparent=true)
close()
