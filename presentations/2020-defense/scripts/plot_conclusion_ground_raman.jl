#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using PyCall
using DataStructures
using MAT

const data = matread(joinpath(@__DIR__, "../data/groundstate_raman_rabi.mat"))

const x = data["x"]
const y = data["y"]
const dy = data["dy"]

function model_rabi(x, p)
    offset, amp, Γ, γ, Ω = p
    return offset .+ amp .* exp.(.-Γ .* x) .* (1 .- (1 .- exp.(.-γ .* x) .* cos.(Ω .* x)) ./ 2)
end

# const param = [0.3833, 0.287, 0.0133, 0.067, 2π * 0.1705]
const param = [0, 1, 0.0133, 0.067, 2π * 0.1705]

const prefix = joinpath(@__DIR__, "../figures/conclusion_ground_raman")

const plotx = linspace(0, 26, 1000)

figure()
errorbar(x[2:end] .* 1e6, (y[2:end] .- 0.3833) ./ 0.287,
         dy[2:end] ./ 0.287, fmt="C0o")
plot(plotx, model_rabi(plotx, param), "C0")
grid()
xlim([0, 27])
ylim([0, 1.0])
xlabel("Raman Time (\$\\mathrm{\\mu s}\$)")
ylabel("Fraction of Feshbach molecule")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
