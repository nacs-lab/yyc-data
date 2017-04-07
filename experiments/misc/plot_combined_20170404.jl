#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const iname_a = joinpath(@__DIR__, "data", "data_20170402_205344.csv")
const params_a, ratios_a, uncs_a = NaCsData.calc_survival(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20170404_133229.csv")
const params_b, ratios_b, uncs_b = NaCsData.calc_survival(iname_b)

# With cooling +-1
Params_A1_1 = [linspace(-18.985, -18.785, 11);] .* 1000 .+ 18462.5
Params_A1_2 = [linspace(-18.15, -17.95, 11);] .* 1000 .+ 18462.5
Params_A2_1 = [linspace(-18.945, -19.145, 11);] .* 1000 .+ 18485
Params_A2_2 = [linspace(-18.00, -17.75, 11);] .* 1000 .+ 18485
Params_A3_1 = [linspace(-18.545, -18.585, 11);] .* 1000 .+ 18496.5
Params_A3_2 = [linspace(-18.415, -18.455, 11);] .* 1000 .+ 18496.5

# Without cooling +-1, -2
Params_A4_1 = [linspace(-18.985, -18.785, 11);] .* 1000 .+ 18462.5
Params_A4_2 = [linspace(-18.15, -17.95, 11);] .* 1000 .+ 18462.5
Params_A4_3 = [linspace(-17.53, -17.745, 11);] .* 1000 .+ 18462.5
Params_A5_1 = [linspace(-18.945, -19.145, 11);] .* 1000 .+ 18485
Params_A5_2 = [linspace(-18.00, -17.75, 11);] .* 1000 .+ 18485
Params_A5_3 = [linspace(-17.14, -17.44, 11);] .* 1000 .+ 18485
Params_A6_1 = [linspace(-18.545, -18.585, 11);] .* 1000 .+ 18496.5
Params_A6_2 = [linspace(-18.415, -18.455, 11);] .* 1000 .+ 18496.5
Params_A6_3 = [linspace(-18.35, -18.39, 11);] .* 1000 .+ 18496.5

# Without cooling carrier
Params_A7 = [linspace(-18.335, -18.58, 11);] .* 1000 .+ 18462.5
Params_A8 = [linspace(-18.31, -18.61, 11);] .* 1000 .+ 18485
Params_A9 = [linspace(-18.47, -18.53, 11);] .* 1000 .+ 18496.5

# Without cooling axial high orders
Params_A10 = [linspace(-18.35, -17.90, 91);] .* 1000 .+ 18496.5

# Without repeating
Params_A11_1 = [linspace(-18.985, -18.785, 11);] .* 1000 .+ 18462.5
Params_A11_2 = [linspace(-18.15, -17.95, 11);] .* 1000 .+ 18462.5
Params_A12_1 = [linspace(-18.945, -19.145, 11);] .* 1000 .+ 18485
Params_A12_2 = [linspace(-18.00, -17.75, 11);] .* 1000 .+ 18485
Params_A13_1 = [linspace(-18.545, -18.585, 11);] .* 1000 .+ 18496.5
Params_A13_2 = [linspace(-18.415, -18.455, 11);] .* 1000 .+ 18496.5

# With waiting
Params_A14_1 = [linspace(-18.985, -18.785, 11);] .* 1000 .+ 18462.5
Params_A14_2 = [linspace(-18.15, -17.95, 11);] .* 1000 .+ 18462.5
Params_A15_1 = [linspace(-18.945, -19.145, 11);] .* 1000 .+ 18485
Params_A15_2 = [linspace(-18.00, -17.75, 11);] .* 1000 .+ 18485
Params_A16_1 = [linspace(-18.545, -18.585, 11);] .* 1000 .+ 18496.5
Params_A16_2 = [linspace(-18.415, -18.455, 11);] .* 1000 .+ 18496.5

Params_B1 = [linspace(-17.53, -17.745, 11);] .* 1000 .+ 18462.5
Params_B2 = [linspace(-17.14, -17.44, 11);] .* 1000 .+ 18496.5
Params_B3 = [linspace(-18.40, -17.90, 51);] .* 1000 .+ 18496.5
Params_B4 = [linspace(-18.40, -18.35, 11);] .* 1000 .+ 18496.5
Params_B5 = [0:2:100;]
Params_B6 = [2:2:100;]
Params_B7 = [10:10:400;]
Params_B8 = [2:2:30;]
Params_B9 = [2:2:30;]
Params_B10 = [10:10:140;]

offset_a1_1 = length(Params_A1_1)
offset_a1_2 = offset_a1_1 + length(Params_A1_2)
offset_a2_1 = offset_a1_2 + length(Params_A2_1)
offset_a2_2 = offset_a2_1 + length(Params_A2_2)
offset_a3_1 = offset_a2_2 + length(Params_A3_1)
offset_a3_2 = offset_a3_1 + length(Params_A3_2)
offset_a4_1 = offset_a3_2 + length(Params_A4_1)
offset_a4_2 = offset_a4_1 + length(Params_A4_2)
offset_a4_3 = offset_a4_2 + length(Params_A4_3)
offset_a5_1 = offset_a4_3 + length(Params_A5_1)
offset_a5_2 = offset_a5_1 + length(Params_A5_2)
offset_a5_3 = offset_a5_2 + length(Params_A5_3)
offset_a6_1 = offset_a5_3 + length(Params_A6_1)
offset_a6_2 = offset_a6_1 + length(Params_A6_2)
offset_a6_3 = offset_a6_2 + length(Params_A6_3)
offset_a7 = offset_a6_3 + length(Params_A7)
offset_a8 = offset_a7 + length(Params_A8)
offset_a9 = offset_a8 + length(Params_A9)
offset_a10 = offset_a9 + length(Params_A10)
offset_a11_1 = offset_a10 + length(Params_A11_1)
offset_a11_2 = offset_a11_1 + length(Params_A11_2)
offset_a12_1 = offset_a11_2 + length(Params_A12_1)
offset_a12_2 = offset_a12_1 + length(Params_A12_2)
offset_a13_1 = offset_a12_2 + length(Params_A13_1)
offset_a13_2 = offset_a13_1 + length(Params_A13_2)
offset_a14_1 = offset_a13_2 + length(Params_A14_1)
offset_a14_2 = offset_a14_1 + length(Params_A14_2)
offset_a15_1 = offset_a14_2 + length(Params_A15_1)
offset_a15_2 = offset_a15_1 + length(Params_A15_2)
offset_a16_1 = offset_a15_2 + length(Params_A16_1)
offset_a16_2 = offset_a16_1 + length(Params_A16_2)
offset_b1 = length(Params_B1)
offset_b2 = offset_b1 + length(Params_B2)
offset_b3 = offset_b2 + length(Params_B3)
offset_b4 = offset_b3 + length(Params_B4)
offset_b5 = offset_b4 + length(Params_B5)
offset_b6 = offset_b5 + length(Params_B6)
offset_b7 = offset_b6 + length(Params_B7)
offset_b8 = offset_b7 + length(Params_B8)
offset_b9 = offset_b8 + length(Params_B9)
offset_b10 = offset_b9 + length(Params_B10)

Idx_A1_1 = [1:offset_a1_1;]
Idx_A1_2 = [(offset_a1_1 + 1):offset_a1_2;]
Idx_A2_1 = [(offset_a1_2 + 1):offset_a2_1;]
Idx_A2_2 = [(offset_a2_1 + 1):offset_a2_2;]
Idx_A3_1 = [(offset_a2_2 + 1):offset_a3_1;]
Idx_A3_2 = [(offset_a3_1 + 1):offset_a3_2;]
Idx_A4_1 = [(offset_a3_2 + 1):offset_a4_1;]
Idx_A4_2 = [(offset_a4_1 + 1):offset_a4_2;]
Idx_A4_3 = [(offset_a4_2 + 1):offset_a4_3;]
Idx_A5_1 = [(offset_a4_3 + 1):offset_a5_1;]
Idx_A5_2 = [(offset_a5_1 + 1):offset_a5_2;]
Idx_A5_3 = [(offset_a5_2 + 1):offset_a5_3;]
Idx_A6_1 = [(offset_a5_3 + 1):offset_a6_1;]
Idx_A6_2 = [(offset_a5_1 + 1):offset_a6_2;]
Idx_A6_3 = [(offset_a5_2 + 1):offset_a6_3;]
Idx_A7 = [(offset_a6_3 + 1):offset_a7;]
Idx_A8 = [(offset_a7 + 1):offset_a8;]
Idx_A9 = [(offset_a8 + 1):offset_a9;]
Idx_A10 = [(offset_a9 + 1):offset_a10;]
Idx_A11_1 = [(offset_a10 + 1):offset_a11_1;]
Idx_A11_2 = [(offset_a11_1 + 1):offset_a11_2;]
Idx_A12_1 = [(offset_a11_2 + 1):offset_a12_1;]
Idx_A12_2 = [(offset_a12_1 + 1):offset_a12_2;]
Idx_A13_1 = [(offset_a12_2 + 1):offset_a13_1;]
Idx_A13_2 = [(offset_a13_1 + 1):offset_a13_2;]
Idx_A14_1 = [(offset_a13_2 + 1):offset_a14_1;]
Idx_A14_2 = [(offset_a14_1 + 1):offset_a14_2;]
Idx_A15_1 = [(offset_a14_2 + 1):offset_a15_1;]
Idx_A15_2 = [(offset_a15_1 + 1):offset_a15_2;]
Idx_A16_1 = [(offset_a15_2 + 1):offset_a16_1;]
Idx_A16_2 = [(offset_a16_1 + 1):offset_a16_2;]
Idx_B1 = [1:offset_b1;]
Idx_B2 = [(offset_b1 + 1):offset_b2;]
Idx_B3 = [(offset_b2 + 1):offset_b3;]
Idx_B4 = [(offset_b3 + 1):offset_b4;]
Idx_B5 = [(offset_b4 + 1):offset_b5;]
Idx_B6 = [(offset_b4 + 1); (offset_b5 + 1):offset_b6]
Idx_B7 = [(offset_b4 + 1); (offset_b6 + 1):offset_b7]
Idx_B8 = [(offset_b4 + 1); (offset_b7 + 1):offset_b8]
Idx_B9 = [(offset_b4 + 1); (offset_b8 + 1):offset_b9]
Idx_B10 = [(offset_b4 + 1); (offset_b9 + 1):offset_b10]

unshift!(Params_B6, 0)
unshift!(Params_B7, 0)
unshift!(Params_B8, 0)
unshift!(Params_B9, 0)
unshift!(Params_B10, 0)

function plot_params(ratios, uncs, Params, Idx, scale=1; kws...)
    perm = sortperm(Params)
    Params = Params[perm]
    Idx = Idx[perm]
    Ratios = ratios[Idx, 2] .* scale
    Uncs = uncs[Idx, 2] .* scale
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

const prefix = joinpath(@__DIR__, "imgs", "data_combined_20170404")

figure()
# Without cooling
plot_params(ratios_a, uncs_a, Params_A4_1, Idx_A4_1, 1 / 0.95, fmt="ro-", label="Before")
plot_params(ratios_a, uncs_a, Params_A4_2, Idx_A4_2, 1 / 0.95, fmt="ro-")
plot_params(ratios_a, uncs_a, Params_A4_3, Idx_A4_3, 1 / 0.95, fmt="ro-")
plot_params(ratios_a, uncs_a, Params_A7, Idx_A7, 1 / 0.95, fmt="ro-")
# With cooling
plot_params(ratios_a, uncs_a, Params_A1_1, Idx_A1_1, 1 / 0.85, fmt="bo-", label="After")
plot_params(ratios_a, uncs_a, Params_A1_2, Idx_A1_2, 1 / 0.85, fmt="bo-")
plot_params(ratios_b, uncs_b, Params_B1, Idx_B1, 1 / 0.85, fmt="bo-")
grid()
ylim([0, 1])
title("Radial 2")
xlabel("Detuning from carrier (kHz)")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r2")

figure()
# Without cooling
plot_params(ratios_a, uncs_a, Params_A5_1, Idx_A5_1, 1 / 0.95, fmt="ro-", label="Before")
plot_params(ratios_a, uncs_a, Params_A5_2, Idx_A5_2, 1 / 0.95, fmt="ro-")
plot_params(ratios_a, uncs_a, Params_A5_3, Idx_A5_3, 1 / 0.95, fmt="ro-")
plot_params(ratios_a, uncs_a, Params_A8, Idx_A8, 1 / 0.95, fmt="ro-")
# With cooling
plot_params(ratios_a, uncs_a, Params_A2_1, Idx_A2_1, 1 / 0.85, fmt="bo-", label="After")
plot_params(ratios_a, uncs_a, Params_A2_2, Idx_A2_2, 1 / 0.85, fmt="bo-")
plot_params(ratios_b, uncs_b, Params_B2, Idx_B2, 1 / 0.85, fmt="bo-")
grid()
ylim([0, 0.9])
title("Radial 3")
xlabel("Detuning from carrier (kHz)")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_r3")

figure(figsize=[2.5, 1] * 4.8)
# Without cooling
plot_params(ratios_a, uncs_a, Params_A6_1, Idx_A6_1, 1 / 0.95, fmt="ro-", label="Before")
plot_params(ratios_a, uncs_a, Params_A6_2, Idx_A6_2, 1 / 0.95, fmt="ro-")
plot_params(ratios_a, uncs_a, Params_A9, Idx_A9, 1 / 0.95, fmt="ro-")
let
    params_plot = [Params_B4; Params_A10[2:end]]
    local idx_a10 = Idx_A10[2:end]
    ratios_plot = [ratios_b[Idx_B4, 2]; ratios_a[idx_a10, 2]] ./ 0.95
    uncs_plot = [uncs_b[Idx_B4, 2]; uncs_a[idx_a10, 2]] ./ 0.95
    errorbar(params_plot, ratios_plot, uncs_plot,
             fmt="r^-", label="Before")
end
# With cooling
plot_params(ratios_a, uncs_a, Params_A3_1, Idx_A3_1, 1 / 0.85, fmt="bo-", label="After")
plot_params(ratios_a, uncs_a, Params_A3_2, Idx_A3_2, 1 / 0.85, fmt="bo-")
plot_params(ratios_b, uncs_b, Params_B3, Idx_B3, 1 / 0.85, fmt="b^-", label="After")
grid()
ylim([0, 0.8])
title("Axial 1")
xlabel("Detuning from carrier (kHz)")
ylabel("Normalized survival")
legend()
maybe_save("$(prefix)_a1")

maybe_show()
