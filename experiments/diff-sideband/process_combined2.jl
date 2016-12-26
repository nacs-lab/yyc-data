#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

# Hard coded for now
const freqs_MHz = linspace(-18.4, -17.9, 51);
# const freqs_MHz = linspace(-18.4, -18.0, 41);
const nfreqs = length(freqs_MHz)

const params_c, ratios_c, uncs_c = NaCsData.calc_survival(ARGS[1])
const params = freqs_MHz
const ratios1_r = ratios_c[3:(nfreqs + 2), 2]
const ratios2_r = ratios_c[(nfreqs + 3):(nfreqs * 2 + 2), 2]
const uncs1_r = uncs_c[3:(nfreqs + 2), 2]
const uncs2_r = uncs_c[(nfreqs + 3):(nfreqs * 2 + 2), 2]

ratios_diff = ratios2_r .- ratios1_r
uncs_diff = sqrt.(uncs2_r.^2 .+ uncs1_r.^2)

factor_cold = ratios_c[1, 2]
unc_factor_cold = uncs_c[1, 2]
factor_hot = ratios_c[2, 2] - ratios_c[1, 2]
unc_factor_hot = sqrt(uncs_c[2, 2]^2 + uncs_c[1, 2]^2)
factor_all = ratios_c[2, 2]
unc_factor_all = uncs_c[2, 2]

ratios_cold = ratios1_r ./ factor_cold
ratios_hot = ratios_diff ./ factor_hot
ratios_all = ratios2_r ./ factor_all
uncs_cold = uncs1_r ./ factor_cold
uncs_hot = uncs_diff ./ factor_hot
uncs_all = uncs2_r ./ factor_all

figure()
errorbar(params, ratios_cold, uncs_cold, label="Cold", fmt="bo-")
errorbar(params, ratios_hot, uncs_hot, label="Hot", fmt="ro-")
legend()
grid()
ylim([0, ylim()[2]])
xlabel("\$f\$ (MHz)")
ylabel("Normalized survival")
savefig("$(ARGS[2])_hot.png", bbox_inches="tight", transparent=true)
close()

figure()
errorbar(params, ratios_all, uncs_all, label="Total", fmt="go-")
errorbar(params, ratios_cold, uncs_cold, label="Cold", fmt="bo-")
legend()
grid()
ylim([0, ylim()[2]])
xlabel("\$f\$ (MHz)")
ylabel("Normalized survival")
savefig("$(ARGS[2])_all.png", bbox_inches="tight", transparent=true)
close()

show()
