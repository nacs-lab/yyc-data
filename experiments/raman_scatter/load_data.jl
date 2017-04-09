#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData

# 1. Co-prop F1: 0.25
# 2. Co-prop F2: 0.22
# 3. Up F1: 1.0
# 4. Down F1: 0.22
# 5. Counter-OP F2: 0.05

const iname_a = joinpath(@__DIR__, "data", "data_20170405_114322.csv")
const params_a, ratios_a, uncs_a = NaCsData.calc_survival(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20170405_151151.csv")
const params_b, ratios_b, uncs_b = NaCsData.calc_survival(iname_b)

Params_A1 = [0, 2, 5, 10, 20, 30, 40, 160]
Params_A2 = [2, 5, 10, 20, 30, 40, 160]
Params_A3 = [2, 5, 10, 20, 40, 80, 160]
Params_A4 = [0.2, 0.5, 1, 2, 5, 7.5, 10, 40]
Params_A5 = [2, 5, 10, 20, 30, 40, 80]
Params_B1 = [240]
Params_B2 = [8]
Params_B3 = [1, 8]
Params_B4 = [80, 160]
Params_B5 = [8, 15]

offset_a1 = length(Params_A1)
offset_a2 = offset_a1 + length(Params_A2)
offset_a3 = offset_a2 + length(Params_A3)
offset_a4 = offset_a3 + length(Params_A4)
offset_a5 = offset_a4 + length(Params_A5)
offset_b1 = length(Params_B1)
offset_b2 = offset_b1 + length(Params_B2)
offset_b3 = offset_b2 + length(Params_B3)
offset_b4 = offset_b3 + length(Params_B4)
offset_b5 = offset_b4 + length(Params_B5)

Idx_A1 = [1:offset_a1;]
Idx_A2 = [1; (offset_a1 + 1):offset_a2]
Idx_A3 = [1; (offset_a2 + 1):offset_a3]
Idx_A4 = [1; (offset_a3 + 1):offset_a4]
Idx_A5 = [1; (offset_a4 + 1):offset_a5]
Idx_B1 = [1:offset_b1;]
Idx_B2 = [(offset_b1 + 1):offset_b2;]
Idx_B3 = [(offset_b2 + 1):offset_b3;]
Idx_B4 = [(offset_b3 + 1):offset_b4;]
Idx_B5 = [(offset_b4 + 1):offset_b5;]

params_f1_coprop = [Params_A1; Params_B1]
params_f2_coprop = [0; Params_A2; Params_B2]
params_f1_up = [0; Params_A3; Params_B3]
params_f1_down = [0; Params_A4; Params_B4]
params_f2_counterop = [0; Params_A5; Params_B5]

ratios_f1_coprop = [ratios_a[Idx_A1, 2]; ratios_b[Idx_B1, 2]]
ratios_f2_coprop = [ratios_a[Idx_A2, 2]; ratios_b[Idx_B2, 2]]
ratios_f1_up = [ratios_a[Idx_A3, 2]; ratios_b[Idx_B3, 2]]
ratios_f1_down = [ratios_a[Idx_A4, 2]; ratios_b[Idx_B4, 2]]
ratios_f2_counterop = [ratios_a[Idx_A5, 2]; ratios_b[Idx_B5, 2]]

uncs_f1_coprop = [uncs_a[Idx_A1, 2]; uncs_b[Idx_B1, 2]]
uncs_f2_coprop = [uncs_a[Idx_A2, 2]; uncs_b[Idx_B2, 2]]
uncs_f1_up = [uncs_a[Idx_A3, 2]; uncs_b[Idx_B3, 2]]
uncs_f1_down = [uncs_a[Idx_A4, 2]; uncs_b[Idx_B4, 2]]
uncs_f2_counterop = [uncs_a[Idx_A5, 2]; uncs_b[Idx_B5, 2]]

function sort_params!(params, ratios, uncs)
    perm = sortperm(params)
    params .= params[perm]
    ratios .= ratios[perm]
    uncs .= uncs[perm]
    return
end

sort_params!(params_f1_coprop, ratios_f1_coprop, uncs_f1_coprop)
sort_params!(params_f2_coprop, ratios_f2_coprop, uncs_f2_coprop)
sort_params!(params_f1_up, ratios_f1_up, uncs_f1_up)
sort_params!(params_f1_down, ratios_f1_down, uncs_f1_down)
sort_params!(params_f2_counterop, ratios_f2_counterop, uncs_f2_counterop)
