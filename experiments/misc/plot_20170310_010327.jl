#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])
with_axial1 = params .<= 0
without_axial1 = params .>= 0

errorbar(params[without_axial1][1:end - 1], ratios[without_axial1, 2][1:end - 1],
         uncs[without_axial1, 2][1:end - 1],
         label="Without axial 1", fmt="ro-")
errorbar(.-params[with_axial1][2:end], ratios[with_axial1, 2][2:end],
         uncs[with_axial1, 2][2:end],
         label="With axial 1", fmt="bo-")
legend()
grid()
xlabel("VTweezer")
ylabel("Survival")
const prefix = ARGS[2]
savefig("$(prefix).png", bbox_inches="tight", transparent=true)
savefig("$(prefix).svg", bbox_inches="tight", transparent=true)
