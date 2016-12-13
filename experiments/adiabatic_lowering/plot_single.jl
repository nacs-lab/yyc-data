#!/usr/bin/julia

include("process.jl")

const prefix = ARGS[1]
const powers = readcsv(ARGS[2])
const params, ratios, uncs = NaCsData.calc_survival(ARGS[3])

((depths, ratios1, uncs1),
 (plot_depth, fit_cdf, fit_pdf)) = process_lowering(powers, params, ratios, uncs)

figure()
errorbar(depths, ratios1, uncs1)
plot(plot_depth, fit_cdf)
ylim([0, ylim()[2]])
xlabel("\$T\$ (MHz)")
ylabel("Survival")
grid()
savefig("$(prefix)_cdf.png", bbox_inches="tight", transparent=true)
close()

figure()
plot(plot_depth, fit_pdf)
ylim([0, ylim()[2]])
xlabel("\$T\$ (MHz)")
ylabel("Population")
grid()
savefig("$(prefix)_pdf.png", bbox_inches="tight", transparent=true)
close()
