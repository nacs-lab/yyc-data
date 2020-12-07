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
using Statistics

using DelimitedFiles

function read_csv_compressed(fname, args...; kwargs...)
    LibArchive.Reader(fname) do reader
        LibArchive.support_format_raw(reader)
        LibArchive.support_filter_all(reader)
        LibArchive.next_header(reader)
        readdlm(reader, args...; kwargs...)
    end
end

const data = read_csv_compressed(
    joinpath(@__DIR__, "../../thesis/data/Timebase1038ASE_20200718_193100/data.csv.zst"),
    ',', Float64)

const params = (
    bw = 2,
    xthresh = 0.00005,
    yzero = 0.006,
    # xcutoff1 = 0.2705,
    # xcutoff2 = 0.301,
    lambda0 = 999.85, # nm
    lambda1 = 1100.85, # nm
)

function process_data(data, params)
    curxs = [data[1, 1]]
    ymin = ymax = data[1, 2]
    xs = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 2:size(data, 1)
        x = data[i, 1]
        y = data[i, 2]
        if x - curxs[end] <= params.xthresh
            push!(curxs, x)
            if y < ymin
                ymin = y
            elseif y > ymax
                ymax = y
            end
            continue
        end
        if length(curxs) > 1
            push!(xs, mean(curxs) * (params.lambda1 - params.lambda0) + params.lambda0)
            if ymin < params.yzero
                ymin = -0.0001
            end
            push!(ys, (ymin + ymax) / 2 / params.bw)
            push!(uncs, (ymax - ymin) / 2 / params.bw)
        end
        resize!(curxs, 1)
        curxs[1] = x
        ymin = ymax = y
    end
    if length(curxs) > 1
        push!(xs, mean(curxs) * (params.lambda1 - params.lambda0) + params.lambda0)
        push!(ys, (ymin + ymax) / 2 / params.bw)
        push!(uncs, (ymax - ymin) / 2 / params.bw)
    end
    return (xs=xs, ys=ys, uncs=uncs)
end

const prefix = joinpath(@__DIR__, "../figures/raman_transfer_spectrum_postfilter")

const processed = process_data(data, params)

figure()
ax = gca()
mask = matplotlib.patches.Rectangle((1037.7, 0), 2.8, 2.4, facecolor="C0", alpha=0.7)
ax.add_patch(mask)
errorbar(processed.xs, processed.ys, processed.uncs, linestyle="None", color="C0")
grid()
xlim([params.lambda0, params.lambda1])
ylim([0, 2.4])
xlabel("Wavelength (nm)")
ylabel("Spectral Density (nW/nm)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
