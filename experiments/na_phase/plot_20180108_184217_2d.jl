#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci
using MAT

const iname_a = joinpath(@__DIR__, "data", "data_20180108_184217.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const counts_p, counts_v, counts_u = matopen(iname_a) do f
    cs = read(f, "Counts")::Array{Float64,3}
    nimgs = size(cs, 1)
    nsites = size(cs, 2)
    data_dict = Dict{Float64,Tuple{Matrix{Float64},Matrix{Float64},typeof(Ref(0))}}()
    for seq in 1:length(params_a)
        param = params_a[seq]
        count = @view cs[:, :, seq]
        if haskey(data_dict, param)
            item = data_dict[param]
        else
            item = data_dict[param] = (zeros(nimgs, nsites), zeros(nimgs, nsites), Ref(0))
        end
        item[1] .+= count
        item[2] .+= count.^2
        item[3][] += 1
    end
    params = sort(collect(keys(data_dict)))
    len = length(params)
    avgs = zeros(nimgs, nsites, len)
    uncs = zeros(nimgs, nsites, len)
    for i in 1:len
        p = data_dict[params[i]]
        nrep = p[3][]
        for j in 1:nsites
            for k in 1:nimgs
                avg = p[1][k, j] / nrep
                avg2 = p[2][k, j] / nrep
                avgs[k, j, i] = avg
                uncs[k, j, i] = sqrt(avg2 - avg^2) / sqrt(nrep - 1)
            end
        end
    end
    params, avgs, uncs
end

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

const nduties = 9
const nphases = 11
const duties_ex = linspace(28, 44, nduties)
const phases_ex = linspace(-30, -130, nphases)
const spec = ((duties_ex for _ in 1:nphases)...)

const data_a = NaCsData.select_count(params_a, logicals_a, selector)
const probs_p, probs_r, probs_u = NaCsData.get_values(data_a)


const split_a = NaCsData.split_data(data_a, spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20180108_184217")

figure()
for i in 1:nphases
    phase = phases_ex[i]
    idx_rng = (nduties * (i - 1) + 1):(nduties * (i - 1) + nduties)
    errorbar(duties_ex, probs_r[idx_rng, 1], probs_u[idx_rng, 1],
             label="$phase\${}^\\circ\$")
end
grid()
ylim([0, 0.6])
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=16)
title("Loading")
xlabel("Duty cycle (%)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_duty_load")

figure()
for i in 1:nphases
    phase = phases_ex[i]
    idx_rng = (nduties * (i - 1) + 1):(nduties * (i - 1) + nduties)
    errorbar(duties_ex, probs_r[idx_rng, 2], probs_u[idx_rng, 2],
             label="$phase\${}^\\circ\$")
end
grid()
ylim([0.5, 1])
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=16)
title("Survival")
xlabel("Duty cycle (%)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_duty_survival")

figure()
for i in 1:nphases
    phase = phases_ex[i]
    idx_rng = (nduties * (i - 1) + 1):(nduties * (i - 1) + nduties)
    errorbar(duties_ex, counts_v[1, 1, idx_rng], counts_u[1, 1, idx_rng],
             label="$phase\${}^\\circ\$")
end
grid()
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=16)
title("Total count first image")
xlabel("Duty cycle (%)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_duty_count1")

figure()
for i in 1:nphases
    phase = phases_ex[i]
    idx_rng = (nduties * (i - 1) + 1):(nduties * (i - 1) + nduties)
    errorbar(duties_ex, counts_v[3, 1, idx_rng], counts_u[3, 1, idx_rng],
             label="$phase\${}^\\circ\$")
end
grid()
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=16)
title("Total count second image")
xlabel("Duty cycle (%)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_duty_count2")

figure()
for i in 1:nduties
    duty = duties_ex[i]
    idx_rng = i:nduties:(nphases * nduties)
    errorbar(phases_ex, probs_r[idx_rng, 1], probs_u[idx_rng, 1],
             label="$duty%")
end
grid()
ylim([0, 0.6])
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=18)
title("Loading")
xlabel("Phase (\${}^\\circ\$)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_phase_load")

figure()
for i in 1:nduties
    duty = duties_ex[i]
    idx_rng = i:nduties:(nphases * nduties)
    errorbar(phases_ex, probs_r[idx_rng, 2], probs_u[idx_rng, 2],
             label="$duty%")
end
grid()
ylim([0.5, 1])
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=18)
title("Survival")
xlabel("Phase (\${}^\\circ\$)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_phase_survival")

figure()
for i in 1:nduties
    duty = duties_ex[i]
    idx_rng = i:nduties:(nphases * nduties)
    errorbar(phases_ex, counts_v[1, 1, idx_rng], counts_u[1, 1, idx_rng],
             label="$duty%")
end
grid()
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=18)
title("Total count first image")
xlabel("Phase (\${}^\\circ\$)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_phase_count1")

figure()
for i in 1:nduties
    duty = duties_ex[i]
    idx_rng = i:nduties:(nphases * nduties)
    errorbar(phases_ex, counts_v[3, 1, idx_rng], counts_u[3, 1, idx_rng],
             label="$duty%")
end
grid()
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=18)
title("Total count second image")
xlabel("Phase (\${}^\\circ\$)")
ylabel("Probability")
NaCsPlot.maybe_save("$(prefix)_phase_count2")

const phases_plot = linspace(-25, -135, 10000)

function post_show_2d()
    colorbar()
    for p0 in -30:30:150
        plot(phases_plot, (44 + p0 * 5 / 9) .+ phases_plot .* (5 / 9),
             color="orange", linewidth=1.5)
        plot(phases_plot, (44 + p0 * 5 / 18) .+ phases_plot .* (5 / 18),
             color="white", linewidth=1.5)
    end
    xlim([-25, -135])
    ylim([45, 27])
    xticks(linspace(-30, -130, 6))
    xlabel("Phase (\${}^\\circ\$)")
    yticks(linspace(28, 44, 5))
    ylabel("Duty cycle (%)")
end

figure()
imshow(reshape(probs_r[:, 1], (nduties, nphases)),
       extent=[-25, -135, 45, 27], aspect="auto", vmin=0, vmax=0.6,
       interpolation="nearest")
title("Loading")
post_show_2d()
NaCsPlot.maybe_save("$(prefix)_2d_load")

figure()
imshow(reshape(probs_r[:, 2], (nduties, nphases)),
       extent=[-25, -135, 45, 27], aspect="auto", vmin=0.5, vmax=1,
       interpolation="nearest")
title("Survival")
post_show_2d()
NaCsPlot.maybe_save("$(prefix)_2d_survival")

figure()
imshow(reshape(counts_v[1, 1, :], (nduties, nphases)),
       extent=[-25, -135, 45, 27], aspect="auto", interpolation="nearest")
title("Total count first image")
post_show_2d()
NaCsPlot.maybe_save("$(prefix)_2d_count1")

figure()
imshow(reshape(counts_v[3, 1, :], (nduties, nphases)),
       extent=[-25, -135, 45, 27], aspect="auto", interpolation="nearest")
title("Total count second image")
post_show_2d()
NaCsPlot.maybe_save("$(prefix)_2d_count2")

NaCsPlot.maybe_show()
