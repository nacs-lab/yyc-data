#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
# using DataStructures
# using LsqFit
using LibArchive

using DelimitedFiles

function read_csv_compressed(fname, args...; kwargs...)
    LibArchive.Reader(fname) do reader
        LibArchive.support_format_raw(reader)
        LibArchive.support_filter_all(reader)
        LibArchive.next_header(reader)
        readdlm(reader, args...; kwargs...)
    end
end

const fname31 = joinpath(@__DIR__, "../../2020-postdoc/3311boundstate.csv.zst")
const fname32 = joinpath(@__DIR__, "../../2020-postdoc/3322boundstate.csv.zst")
const fname42 = joinpath(@__DIR__, "../../2020-postdoc/4422boundstate.csv.zst")
const data31 = read_csv_compressed(fname31, ',', Float64)
const data32 = read_csv_compressed(fname32, ',', Float64)
const data42 = read_csv_compressed(fname42, ',', Float64)

const prefix = joinpath(@__DIR__, "../figures/interaction_shift_atomic_wavefunction")


figure()
plot(data42[:, 1], data42[:, 2] .* 2e6, "C2", label="Na(2,2) Cs(4,4)")
plot(data32[:, 1], data32[:, 2] .* 2e6, "C1", label="Na(2,2) Cs(3,3)")
plot(data31[:, 1], data31[:, 2] .* 2e6, "C3--", label="Na(1,1) Cs(3,3)")
axvline(100, color="C0", ls="dashed", linewidth=3)
grid()
xscale("log")
yscale("log")
xlim([2.7, 3000])
ylim([0.0015, 1450])
xticks([10, 100, 1000], ["10", "100", "1000"])
xlabel("R (Angstroms)", fontsize=20)
ylabel("\$\\left|\\psi(R)\\right|^2\$", fontsize=20)
gca().tick_params(axis="both", which="major", labelsize=15)
legend(fontsize="x-small")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
