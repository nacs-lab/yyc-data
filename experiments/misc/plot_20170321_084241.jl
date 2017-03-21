#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])
without_cool = params .<= 0
with_cool = params .>= 0

errorbar(sqrt.(.-params[without_cool]), ratios[without_cool, 2], uncs[without_cool, 2],
         label="Before", fmt="ro-")
errorbar(sqrt.(params[with_cool]), ratios[with_cool, 2], uncs[with_cool, 2],
         label="After", fmt="bo-")
ylim([0, 1])
legend()
grid()
xlabel("\$\\sqrt{VTweezer}\$")
ylabel("Survival")

const prefix = ARGS[2]
savefig("$(prefix).png", bbox_inches="tight", transparent=true)
savefig("$(prefix).svg", bbox_inches="tight", transparent=true)

close()
# show()
