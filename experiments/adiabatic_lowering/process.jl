#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using Dierckx
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const prefix = ARGS[1]

const powers = readcsv(ARGS[2])
const params, ratios, uncs = NaCsData.calc_survival(ARGS[3])

unc2w(x) = if x <= 0
    return 1.0
else
    return 1 / x
end

function process_lowering(powers, params, ratios, uncs)
    orig = powers[1, 2]
    amps = powers[end:-1:2, 1]
    lower_ratios = powers[end:-1:2, 2] ./ orig
    full_depth = 2 * 20 # in MHz. Hard coded for now
    fit_depths = full_depth .* sqrt.(lower_ratios)
    spl = Spline1D(amps, fit_depths, s=0.1)
    depths = spl(params)
    uncs_1d = uncs[:, 2]
    fit = Spline1D(depths, ratios[:, 2],
                   k=3, w=unc2w.(uncs_1d), s=length(uncs_1d))
    plot_depth = [linspace(0, depths[end] * 0.85, 1000);]
    return ((depths, ratios[:, 2], uncs[:, 2]),
            (plot_depth, fit(plot_depth), derivative(fit, plot_depth)))
end

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
