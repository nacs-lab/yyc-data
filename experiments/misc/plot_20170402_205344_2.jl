#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const iname = joinpath(@__DIR__, "data", "data_20170402_205344.csv")
const params, ratios, uncs = NaCsData.calc_survival(iname)

# With cooling +-1
Params1_1 = [linspace(-18.985, -18.785, 11) + 18.4625;] .* 1000
Params1_2 = [linspace(-18.15, -17.95, 11) + 18.4625;] .* 1000
Params2_1 = [linspace(-18.945, -19.145, 11) + 18.485;] .* 1000
Params2_2 = [linspace(-18.00, -17.75, 11) + 18.485;] .* 1000
Params3_1 = [linspace(-18.545, -18.585, 11) + 18.4965;] .* 1000
Params3_2 = [linspace(-18.415, -18.455, 11) + 18.4965;] .* 1000

# Without cooling +-1, -2
Params4_1 = [linspace(-18.985, -18.785, 11);]
Params4_2 = [linspace(-18.15, -17.95, 11);]
Params4_3 = [linspace(-17.53, -17.745, 11);]
Params5_1 = [linspace(-18.945, -19.145, 11);]
Params5_2 = [linspace(-18.00, -17.75, 11);]
Params5_3 = [linspace(-17.14, -17.44, 11);]
Params6_1 = [linspace(-18.545, -18.585, 11);]
Params6_2 = [linspace(-18.415, -18.455, 11);]
Params6_3 = [linspace(-18.35, -18.39, 11);]

# Without cooling carrier
Params7 = [linspace(-18.335, -18.58, 11);]
Params8 = [linspace(-18.31, -18.61, 11);]
Params9 = [linspace(-18.47, -18.53, 11);]

# Without cooling axial high orders
Params10 = [linspace(-18.35, -17.90, 91);]

# Without repeating
Params11_1 = [linspace(-18.985, -18.785, 11);]
Params11_2 = [linspace(-18.15, -17.95, 11);]
Params12_1 = [linspace(-18.945, -19.145, 11);]
Params12_2 = [linspace(-18.00, -17.75, 11);]
Params13_1 = [linspace(-18.545, -18.585, 11);]
Params13_2 = [linspace(-18.415, -18.455, 11);]

# With waiting
Params14_1 = [linspace(-18.985, -18.785, 11);]
Params14_2 = [linspace(-18.15, -17.95, 11);]
Params15_1 = [linspace(-18.945, -19.145, 11);]
Params15_2 = [linspace(-18.00, -17.75, 11);]
Params16_1 = [linspace(-18.545, -18.585, 11);]
Params16_2 = [linspace(-18.415, -18.455, 11);]

offset1_1 = length(Params1_1)
offset1_2 = offset1_1 + length(Params1_2)
offset2_1 = offset1_2 + length(Params2_1)
offset2_2 = offset2_1 + length(Params2_2)
offset3_1 = offset2_2 + length(Params3_1)
offset3_2 = offset3_1 + length(Params3_2)
offset4_1 = offset3_2 + length(Params4_1)
offset4_2 = offset4_1 + length(Params4_2)
offset4_3 = offset4_2 + length(Params4_3)
offset5_1 = offset4_3 + length(Params5_1)
offset5_2 = offset5_1 + length(Params5_2)
offset5_3 = offset5_2 + length(Params5_3)
offset6_1 = offset5_3 + length(Params6_1)
offset6_2 = offset6_1 + length(Params6_2)
offset6_3 = offset6_2 + length(Params6_3)
offset7 = offset6_3 + length(Params7)
offset8 = offset7 + length(Params8)
offset9 = offset8 + length(Params9)
offset10 = offset9 + length(Params10)
offset11_1 = offset10 + length(Params11_1)
offset11_2 = offset11_1 + length(Params11_2)
offset12_1 = offset11_2 + length(Params12_1)
offset12_2 = offset12_1 + length(Params12_2)
offset13_1 = offset12_2 + length(Params13_1)
offset13_2 = offset13_1 + length(Params13_2)
offset14_1 = offset13_2 + length(Params14_1)
offset14_2 = offset14_1 + length(Params14_2)
offset15_1 = offset14_2 + length(Params15_1)
offset15_2 = offset15_1 + length(Params15_2)
offset16_1 = offset15_2 + length(Params16_1)
offset16_2 = offset16_1 + length(Params16_2)

Idx1_1 = [1:offset1_1;]
Idx1_2 = [(offset1_1 + 1):offset1_2;]
Idx2_1 = [(offset1_2 + 1):offset2_1;]
Idx2_2 = [(offset2_1 + 1):offset2_2;]
Idx3_1 = [(offset2_2 + 1):offset3_1;]
Idx3_2 = [(offset3_1 + 1):offset3_2;]
Idx4_1 = [(offset3_2 + 1):offset4_1;]
Idx4_2 = [(offset4_1 + 1):offset4_2;]
Idx4_3 = [(offset4_2 + 1):offset4_3;]
Idx5_1 = [(offset4_3 + 1):offset5_1;]
Idx5_2 = [(offset5_1 + 1):offset5_2;]
Idx5_3 = [(offset5_2 + 1):offset5_3;]
Idx6_1 = [(offset5_3 + 1):offset6_1;]
Idx6_2 = [(offset5_1 + 1):offset6_2;]
Idx6_3 = [(offset5_2 + 1):offset6_3;]
Idx7 = [(offset6_3 + 1):offset7;]
Idx8 = [(offset7 + 1):offset8;]
Idx9 = [(offset8 + 1):offset9;]
Idx10 = [(offset9 + 1):offset10;]
Idx11_1 = [(offset10 + 1):offset11_1;]
Idx11_2 = [(offset11_1 + 1):offset11_2;]
Idx12_1 = [(offset11_2 + 1):offset12_1;]
Idx12_2 = [(offset12_1 + 1):offset12_2;]
Idx13_1 = [(offset12_2 + 1):offset13_1;]
Idx13_2 = [(offset13_1 + 1):offset13_2;]
Idx14_1 = [(offset13_2 + 1):offset14_1;]
Idx14_2 = [(offset14_1 + 1):offset14_2;]
Idx15_1 = [(offset14_2 + 1):offset15_1;]
Idx15_2 = [(offset15_1 + 1):offset15_2;]
Idx16_1 = [(offset15_2 + 1):offset16_1;]
Idx16_2 = [(offset16_1 + 1):offset16_2;]

function plot_params(Params, Idx; kws...)
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

const prefix = joinpath(@__DIR__, "imgs", "data_20170402_205344")

figure()
# Axial
plot_params(Params3_1, Idx3_1, fmt="bo-")
plot_params(Params3_2, Idx3_2, fmt="bo-")
# Radial 2
plot_params(Params1_1, Idx1_1, fmt="ro-")
plot_params(Params1_2, Idx1_2, fmt="ro-")
# Radial 3
plot_params(Params2_1, Idx2_1, fmt="ko-")
plot_params(Params2_2, Idx2_2, fmt="ko-")
grid()
ylim([0, 0.9])
xlabel("\$\\Delta_{Raman}\$/MHz")
ylabel("F1 population")
maybe_save("$(prefix)_cold-3d")

maybe_show()
