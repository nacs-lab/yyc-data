#!/usr/bin/julia -f

# β, δ, t, <final energy ± unc>, <escape time ± unc>, <photon count ± unc>
data = readcsv("gamma_5_14.csv", Float32)

using PyPlot

# const cmap1 = PyPlot.ColorMap("autumn")
# const cmap2 = PyPlot.ColorMap("winter")
const plot_β1 = 1f0
const plot_β2 = 0.6f0
const cmap1 = PyPlot.ColorMap("autumn")
const cmap2 = PyPlot.ColorMap("cool_r")

function plot_params(name, all_data, val_idx)
    data_d1 = Dict{Float32,NTuple{3,Vector{Float32}}}()
    data_d2 = Dict{Float32,NTuple{3,Vector{Float32}}}()
    for i in 1:size(all_data, 1)
        β = all_data[i, 1]
        δ = all_data[i, 2]
        t = all_data[i, 3]
        # The t=20000 isn't too different from t=10000 and just add more noise
        # to the plot
        val = all_data[i, val_idx]
        unc = all_data[i, val_idx + 1]

        data_dict = if β == plot_β1
            data_d1
        elseif β == plot_β2
            data_d2
        else
            continue
        end

        data = if !(t in keys(data_dict))
            data_dict[t] = (Float32[], Float32[], Float32[])
        else
            data_dict[t]
        end
        push!(data[1], δ)
        push!(data[2], val)
        push!(data[3], unc)
    end
    for (β, data_dict, cmap) in ((plot_β1, data_d1, cmap1),
                                  (plot_β2, data_d2, cmap2))
        ts = collect(keys(data_dict))
        sort!(ts)
        tmin = ts[1]
        tmax = ts[end]
        for t in ts
            color = cmap((log(t) - log(tmin)) / (log(tmax) - log(tmin)))
            data = data_dict[t]
            errorbar(data[1], data[2], data[3], color=color,
                     label="\$\\beta=$β; trapf=$t\$")
        end
    end
    ylabel(name)
    xlabel("Detuning (MHz)")
    legend(loc=3)
    grid()
end

figure()
plot_params("Final Energy", data, 4)
figure()
plot_params("Escape Time", data, 6)
figure()
plot_params("Photon Count", data, 8)
show()
# savefig("gamma_5_beta1-2.5.png")
