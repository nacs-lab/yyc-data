#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

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

const fname31 = joinpath(@__DIR__, "3311boundstate.csv.zst")
const fname32 = joinpath(@__DIR__, "3322boundstate.csv.zst")
const fname42 = joinpath(@__DIR__, "4422boundstate.csv.zst")
const data31 = read_csv_compressed(fname31, ',', Float64)
const data32 = read_csv_compressed(fname32, ',', Float64)
const data42 = read_csv_compressed(fname42, ',', Float64)

const prefix = joinpath(@__DIR__, "imgs", "wavefunction")

figure()
plot(data42[:, 1], data42[:, 2] .* 1e4, "C2", label="Na(2,2) Cs(4,4)")
plot(data32[:, 1], data32[:, 2] .* 1e4, "C1", label="Na(2,2) Cs(3,3)")
plot(data31[:, 1], data31[:, 2] .* 1e4, "C3--", label="Na(1,1) Cs(3,3)")
grid()
xlim([0, 3080])
ylim([0, 7.3])
ylabel("Probability density")
xlabel("R (Angstroms)")
legend(fontsize="x-small")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
