#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const name20170330 = joinpath(@__DIR__, "data/data_20170330_012226.csv")
const params, ratios, uncs = NaCsData.calc_survival(name20170330)

Params1p = linspace(-18.985, -18.785, 11) .+ 18.47
Params1m = linspace(-18.15, -17.95, 11) .+ 18.47
Params2p = linspace(-18.945, -19.145, 11) .+ 18.4725
Params2m = linspace(-18.00, -17.75, 11) .+ 18.4725
Params3p = linspace(-18.55, -18.58, 16) .+ 18.5015
Params3m = linspace(-18.42, -18.45, 16) .+ 18.5015

offset1p = length(Params1p)
offset1m = offset1p + length(Params1m)
offset2p = offset1m + length(Params2p)
offset2m = offset2p + length(Params2m)
offset3p = offset2m + length(Params3p)
offset3m = offset3p + length(Params3m)

Idx1p = [1:offset1p;]
Idx1m = [(offset1p + 1):offset1m;]
Idx2p = [(offset1m + 1):offset2p;]
Idx2m = [(offset2p + 1):offset2m;]
Idx3p = [(offset2m + 1):offset3p;]
Idx3m = [(offset3p + 1):offset3m;]

const name20170319 = joinpath(@__DIR__, "data/data_20170319_213303.csv")
const params_h, ratios_h, uncs_h = NaCsData.calc_survival(name20170319)

Params1_h = [linspace(-18.93, -18.71, 12); linspace(-18.125, -17.425, 29)]
Params2_h = [linspace(-18.87, -19.07, 11); linspace(-17.925, -17.1, 34)]
Params3p_h = linspace(-18.93, -18.71, 12) .+ 18.3925
Params3m_h = linspace(-18.125, -17.425, 29) .+ 18.3925
Params4p_h = linspace(-18.87, -19.07, 11) .+ 18.3875
Params4m_h = linspace(-17.925, -17.1, 34) .+ 18.3875

offset1_h = length(Params1_h)
offset2_h = offset1_h + length(Params2_h)
offset3p_h = offset2_h + length(Params3p_h)
offset3m_h = offset3p_h + length(Params3m_h)
offset4p_h = offset3m_h + length(Params4p_h)
offset4m_h = offset4p_h + length(Params4m_h)

Idx1_h = [1:offset1_h;]
Idx2_h = [(offset1_h + 1):offset2_h;]
Idx3p_h = [(offset2_h + 1):offset3p_h;]
Idx3m_h = [(offset3p_h + 1):offset3m_h;]
Idx4p_h = [(offset3m_h + 1):offset4p_h;]
Idx4m_h = [(offset4p_h + 1):offset4m_h;]

const name20170320 = joinpath(@__DIR__, "data/data_20170320_130307-155238.csv")
const params_ha, ratios_ha, uncs_ha = NaCsData.calc_survival(name20170320)

function plot_params(ratios, uncs, Params, Idx; kws...)
    perm = sortperm(Params)
    Params = Params[perm]
    Idx = Idx[perm]
    Ratios = ratios[Idx, 2]
    Uncs = uncs[Idx, 2]
    errorbar(Params, Ratios, Uncs; kws...)
end

const save_fig = get(ENV, "NACS_SAVE_FIG", "true") == "true"

function maybe_save(name)
    if save_fig
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if !save_fig
        show()
    end
end

const prefix = joinpath(@__DIR__, "imgs/compare_20170330")

figure()
plot_params(ratios_h, uncs_h, Params3p_h, Idx3p_h, fmt="ro-", label="Before")
plot_params(ratios_h, uncs_h, Params3m_h, Idx3m_h, fmt="ro-")
plot_params(ratios, uncs, Params1p, Idx1p, fmt="bo-", label="After")
plot_params(ratios, uncs, Params1m, Idx1m, fmt="bo-")
grid()
ylim([0, 0.85])
title("Axis 2 (Radial)")
xlabel("\$\\delta\$/MHz")
ylabel("Survival")
legend()
maybe_save("$(prefix)_r2")

figure()
plot_params(ratios_h, uncs_h, Params4p_h, Idx4p_h, fmt="ro-", label="Before")
plot_params(ratios_h, uncs_h, Params4m_h, Idx4m_h, fmt="ro-")
plot_params(ratios, uncs, Params2p, Idx2p, fmt="bo-", label="After")
plot_params(ratios, uncs, Params2m, Idx2m, fmt="bo-")
grid()
ylim([0, 0.8])
title("Axis 3 (Radial)")
xlabel("\$\\delta\$/MHz")
ylabel("Survival")
legend()
maybe_save("$(prefix)_r3")

figure()
errorbar(params_ha .* 1e-6 .+ 18.4975, ratios_ha[:, 2], uncs_ha[:, 2], fmt="ro-", label="Before")
plot_params(ratios, uncs, Params3p, Idx3p, fmt="bo-", label="After")
plot_params(ratios, uncs, Params3m, Idx3m, fmt="bo-")
grid()
ylim([0, ylim()[2]])
title("Axis 1 (Axial)")
xlabel("\$\\delta\$/MHz")
ylabel("Survival")
legend()
maybe_save("$(prefix)_a1")

maybe_show()
