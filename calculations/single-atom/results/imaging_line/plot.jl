#!/usr/bin/julia -f

# β, δ, t, <final energy ± unc>, <escape time ± unc>, <photon count ± unc>
data = readcsv("gamma_5.csv")

using PyPlot

# const cmap1 = PyPlot.ColorMap("autumn")
# const cmap2 = PyPlot.ColorMap("winter")
const cmap1 = PyPlot.ColorMap("autumn_r")
const cmap2 = PyPlot.ColorMap("cool")

function plot_params(name, all_data, val_idx)
    data06 = Dict{Float32,NTuple{3,Vector{Float32}}}()
    data10 = Dict{Float32,NTuple{3,Vector{Float32}}}()
    for i in 1:size(all_data, 1)
        β = all_data[i, 1]
        δ = all_data[i, 2]
        t = all_data[i, 3]
        # The t=20000 isn't too different from t=10000 and just add more noise
        # to the plot
        t == 20000 && continue
        val = all_data[i, val_idx]
        unc = all_data[i, val_idx + 1]

        data_dict = β == 1 ? data10 : data06
        data = if !(t in keys(data_dict))
            data_dict[t] = (Float32[], Float32[], Float32[])
        else
            data_dict[t]
        end
        push!(data[1], δ)
        push!(data[2], val)
        push!(data[3], unc)
    end
    for (β, data_dict, cmap) in ((0.6, data06, cmap1), (1.0, data10, cmap2))
        ts = collect(keys(data_dict))
        sort!(ts)
        tmin = ts[1]
        tmax = ts[end]
        for t in ts
            color = cmap((log(t) - log(tmin)) / (log(tmax) - log(tmin)))
            data = data_dict[t]
            errorbar(data[1], data[2], data[3], color=color,
                     label="\$\\beta=$β; t=$(round(Int, t))\$")
        end
    end
    ylabel(name)
    xlabel("Detuning (MHz)")
    legend()
    grid()
end

# plot_params("Final Energy", data, 4)
# plot_params("Escape Time", data, 6)
plot_params("Photon Count", data, 8)
savefig("gamma_5.png")