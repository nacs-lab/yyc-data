#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using DelimitedFiles

const iname_pos = joinpath(@__DIR__, "data", "raman_agl_pos_20181124.csv")
const rdata_pos = readdlm(iname_pos, ',', Float64, skipstart=1)
const iname_neg = joinpath(@__DIR__, "data", "raman_agl_neg_20181124.csv")
const rdata_neg = readdlm(iname_neg, ',', Float64, skipstart=1)
const rdata = [rdata_neg; rdata_pos]

const agls = rdata[:, 1]
const res_freqs = (rdata[:, 2] .- 298) .* 1000
const res_freqs_unc = rdata[:, 3] .* 1000

function model(x, p)
    return p[1] .- p[2] .* sind.(x .- p[3]).^2
end

fit = curve_fit(model, agls, res_freqs, res_freqs_unc.^-(2/3),
                [204.0, 35.0, 0.0])
uncs = Unc.(fit.param, estimate_errors(fit))
@show uncs

const prefix = joinpath(@__DIR__, "imgs", "fit_raman_agl_20181124")

plotagl = linspace(-100, 100, 10001)

figure()
plot(plotagl, model.(plotagl, (fit.param,)), "C1")
errorbar(agls, res_freqs, res_freqs_unc, fmt="C0.")
text(-50, 172, "\$f_0=$(uncs[1])\\ \$kHz\n" *
     "\$\\delta=$(uncs[2])\\ \$kHz\n" *
     "\$\\theta_0=$(uncs[3])\\ ^\\circ\$", fontsize="small")
text(-60, 207, "\$f=f_0-\\delta\\sin^2(\\theta-\\theta_0)\$")
grid()
xlim([-110, 110])
ylim([168, 212])
title("Angle dependent light shift")
xlabel("B field angle (\${}^\\circ\$)")
ylabel("Raman resonance (298XXX kHz)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
