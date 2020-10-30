#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

# import NaCsCalc.Format: Unc, Sci
# using NaCsCalc
# using NaCsCalc.Utils: interactive
# using NaCsData
# using NaCsData.Fitting: fit_data, fit_survival
# using NaCsPlot
# using PyPlot
# using DataStructures
# using LsqFit

using LibArchive
using DelimitedFiles

const fname = joinpath(@__DIR__, "data/raw.csv.zst")
const data = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end

# Assignments of the Na/Cs spin states and their energies
# Na, Cs
# (0, 0)   0
# (1, 0)  20
# (0, 1)  24
# (2, 0)  40
# (1, 1)  44
# (0, 2)  48
# (3, 0)  60
# (2, 1)  64
# (1, 2)  68
# (0, 3)  72
# (4, 0)  80
# (3, 1)  84
# (2, 2)  88
# (1, 3)  92
# (0, 4)  96
# (5, 0) 100
# (4, 1) 104
# (3, 2) 108
# (2, 3) 112
# (1, 4) 116

# The index for a=0, we know the parity and assignment of this state
const idx_a0 = first(searchsorted(data[:, 1], 0))
const parity0 = (0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)

function propagate_state(data, idx0, direction)
    δi = direction ? 1 : -1
    lastidx = direction ? size(data, 1) : 1
    nstates = size(data, 2) - 1
    as = data[:, 1]
    cur = data[idx0, 2:end]
    prev = data[idx0 - δi, 2:end]
    # state_reindex[si] is the index of the corresponding state at a=0
    state_reindex = [1:nstates;]
    for idx in idx0 + δi:δi:lastidx
        prediction = cur .+ (cur .- prev) .* ((as[idx] - as[idx - δi]) /
                                              (as[idx - δi] - as[idx - 2δi]))
        si = 1
        while si < nstates
            si += 1
            # See if there's a potential crossing between states
            # by checking if swaping the index leads to straighter lines
            prediction1 = prediction[state_reindex[si - 1]]
            prediction2 = prediction[state_reindex[si]]
            err0 = (abs(prediction1 - data[idx, si]) +
                    abs(prediction2 - data[idx, si + 1]))
            err1 = (abs(prediction2 - data[idx, si]) +
                    abs(prediction1 - data[idx, si + 1]))
            if err1 >= err0
                continue
            end
            # If there is potential crossing,
            # check if the two states have different parity and therefore are allowed to cross
            parity1 = parity0[state_reindex[si - 1]]
            parity2 = parity0[state_reindex[si]]
            if parity1 == parity2
                continue
            end
            # println("swap $si @$idx")
            # Swap the states
            state_reindex[si - 1], state_reindex[si] = state_reindex[si], state_reindex[si - 1]
            si += 1
        end
        data[idx, state_reindex .+ 1] .= data[idx, 2:end]
        prev = cur
        cur = data[idx, 2:end]
    end
end
propagate_state(data, idx_a0, true)
propagate_state(data, idx_a0, false)
# for i in 1:20
#     if parity0[i] == 1
#         plot(data[:, 1], data[:, 1 + i], ls="dotted")
#     else
#         plot(data[:, 1], data[:, 1 + i])
#     end
# end
# grid()
# show()
writedlm(joinpath(@__DIR__, "data/sorted.csv"), data, ',')
