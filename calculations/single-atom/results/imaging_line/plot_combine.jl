#!/usr/bin/julia -f

# β, δ, t, <final energy ± unc>, <escape time ± unc>, <photon count ± unc>
# const global_t = 10000
const global_t = 8000
data = readcsv("gamma_5_$(global_t)_comb.csv")

using PyPlot

const cmap = PyPlot.ColorMap("brg_r")

function plot_params(name, all_data, val_idx)
    # {β: {δs: (vals, uncs)}}
    data_dict = Dict{Float64,Dict{Float64,NTuple{2,Float64}}}()
    for i in 1:size(all_data, 1)
        β = all_data[i, 1]
        δ = all_data[i, 2]
        t = all_data[i, 3]

        val = all_data[i, val_idx]
        unc = all_data[i, val_idx + 1]
        data = if !(β in keys(data_dict))
            data_dict[β] = Dict{Float64,NTuple{2,Float64}}()
        else
            data_dict[β]
        end
        if !(δ in keys(data))
            data[δ] = (val, unc)
        else
            oldval, oldunc = data[δ]
            oldw = 1 / oldunc^2
            w = 1 / unc^2
            newval = (oldval * oldw + val * w) / (oldw + w)
            newunc = 1 / sqrt(oldw + w)
            data[δ] = (newval, newunc)
        end
    end
    βs = collect(keys(data_dict))
    sort!(βs)
    βmin = βs[1]
    βmax = βs[end]
    for β in βs
        color = cmap((β - βmin) / (βmax - βmin))
        data = data_dict[β]
        δs = collect(keys(data))
        sort!(δs)
        vals = similar(δs)
        uncs = similar(δs)
        for i in 1:length(δs)
            vals[i], uncs[i] = data[δs[i]]
        end
        errorbar(δs, vals, uncs, color=color, label="\$\\beta=$β\$")
    end
    title("\$t=$(global_t)\$")
    ylabel(name)
    xlabel("Detuning (MHz)")
    # legend()
    legend(loc=2)
    grid()
end

# figure()
# plot_params("Final Energy", data, 4)
# figure()
# plot_params("Escape Time", data, 6)
figure()
plot_params("Photon Count", data, 8)
# show()
savefig("gamma_5_comb_$(global_t).png")
