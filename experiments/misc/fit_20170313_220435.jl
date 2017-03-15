#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
import NaCsCalc.Format: Unc, Sci
using PyPlot
using LsqFit
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

Params1 = [0:40:400; 1600:40:2000]
offset1 = length(Params1)
Params2 = [linspace(-18.03, -18.08, 11); linspace(-18.942, -18.892, 11)]
offset2 = offset1 + length(Params2)
Params3 = [20:20:200; 800:20:1000]
offset3 = offset2 + length(Params3)
Params4 = [linspace(-17.855, -17.905, 11); linspace(-19.073, -19.123, 11)]
offset4 = offset3 + length(Params4)
Params5 = [20:20:200; 800:20:1000]
offset5 = offset4 + length(Params5)
Params6_1 = linspace(-18.54, -18.57, 16)
offset6_1 = offset5 + length(Params6_1)
Params6_2 = linspace(-18.40, -18.43, 16)
offset6_2 = offset6_1 + length(Params6_2)

Idx1 = [1:offset1;]
Idx2 = [(offset1 + 1):offset2;]
Idx3 = [1; (offset2 + 1):offset3]
Idx4 = [(offset3 + 1):offset4;]
Idx5 = [1; (offset4 + 1):offset5]
Idx6_1 = [(offset5 + 1):offset6_1;]
Idx6_2 = [(offset6_1 + 1):offset6_2;]
unshift!(Params3, 0)
unshift!(Params5, 0)

function model(x, p)
    return p[1] .* exp.(-((x .- p[2]) ./ p[3]).^2)
end

function fit_params(Params, Idx; kws...)
    perm = sortperm(Params)
    Params = Params[perm]
    Idx = Idx[perm]
    Ratios = ratios[Idx, 2]
    Uncs = uncs[Idx, 2]
    errorbar(Params, Ratios, Uncs; kws...)
    minparams = minimum(Params)
    maxparams = maximum(Params)
    xplot = linspace(minparams, maxparams, 1000)
    fit = curve_fit(model, Params, Ratios, [maximum(Ratios), mean(Params), std(Params)])
    plot(xplot, model.(xplot, (fit.param,)))
    h, center, width = Unc.(fit.param, estimate_errors(fit), Sci) .* (1, 1, 1000)
    h, center, width
end

const prefix = ARGS[2]

figure()
height_h, center_h, width_h = fit_params(Params6_1, Idx6_1, fmt="o")
height_c, center_c, width_c = fit_params(Params6_2, Idx6_2, fmt="o")
text(-18.57, 0.27, "Heating:\n  Height=\$$height_h\$\n  Center=\$$center_h\$MHz\n  Width=\$$width_h\$kHz",
     size=15)
text(-18.50, 0.10, "Cooling:\n  Height=\$$height_c\$\n  Center=\$$center_c\$MHz\n  Width=\$$width_c\$kHz",
     size=15)
grid()
ylim([0, 0.4])
title("Axial 1 (Gaussian fit)")
xlabel("f/\$MHz\$")
savefig("$(prefix)_a1_fit.png", bbox_inches="tight", transparent=true)
savefig("$(prefix)_a1_fit.svg", bbox_inches="tight", transparent=true)
close()
