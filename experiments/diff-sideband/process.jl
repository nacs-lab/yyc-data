#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])
const params2, ratios2, uncs2 = NaCsData.calc_survival(ARGS[2])
@assert params == params2

ratios_diff = ratios .- ratios2
uncs_diff = sqrt.(uncs.^2 .+ uncs2.^2)

# errorbar(params, ratios[:, 2], uncs[:, 2], label="Full")
errorbar(params, ratios2[:, 2], uncs2[:, 2], label="Cold", fmt="b")
errorbar(params, ratios_diff[:, 2], uncs_diff[:, 2], label="Hot", fmt="r")
legend()
grid()
ylim([0, ylim()[2]])
savefig(ARGS[3])
# show()
