#!/usr/bin/julia

using DelimitedFiles

const data = readdlm(ARGS[1], ',', Float64)

const PAAmp = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
const PADPAmp = 0:0.05:1

function process(data)
    nparams = length(PAAmp) * length(PADPAmp)
    dt = 0.5e-3
    dead_t = 2e-6
    counts = zeros(Int, nparams)
    sums = zeros(nparams)
    mins = fill(Inf, nparams)
    maxs = fill(-Inf, nparams)
    for i in 1:size(data, 1)
        t = data[i, 1]
        idx = floor(Int, t / dt)
        if idx < 0
            continue
        elseif idx >= nparams
            break
        end
        t -= idx * dt
        if t <= dead_t || t >= dt - dead_t
            continue
        end
        idx += 1
        v = data[i, 3]
        counts[idx] += 1
        sums[idx] += v
        mins[idx] = min(mins[idx], v)
        maxs[idx] = max(maxs[idx], v)
    end
    res = zeros(nparams, 3)
    for i in 1:length(PADPAmp)
        dp = PADPAmp[i]
        for j in 1:length(PAAmp)
            pa = PAAmp[j]
            k = (i - 1) * length(PAAmp) + j
            c = counts[k]
            @assert(c > 6000)
            res[k, 1] = pa
            res[k, 2] = dp
            res[k, 3] = sums[k] / c
            @assert(maxs[k] - mins[k] <= 0.48)
        end
    end
    return res
end

open(ARGS[2], "w") do fh
    write(fh, "PAAOM/AMP,PADPaom/AMP,Power\n")
    writedlm(fh, process(data), ',')
end
