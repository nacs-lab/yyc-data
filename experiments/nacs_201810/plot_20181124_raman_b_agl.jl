#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181124_093048.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))
const iname_b = joinpath(@__DIR__, "data", "data_20181124_104133.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
data_nacs_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3, 4)))
const iname_c = joinpath(@__DIR__, "data", "data_20181124_105909.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
data_nacs_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (3, 4)))
const iname_d = joinpath(@__DIR__, "data", "data_20181124_112626.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)
data_nacs_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, 2), (3, 4)))
const iname_e = joinpath(@__DIR__, "data", "data_20181124_120250.mat")
const params_e, logicals_e = NaCsData.load_striped_mat(iname_e)
data_nacs_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, 2), (3, 4)))
const iname_f = joinpath(@__DIR__, "data", "data_20181124_122205.mat")
const params_f, logicals_f = NaCsData.load_striped_mat(iname_f)
data_nacs_f = NaCsData.select_count(params_f, logicals_f, NaCsData.select_single((1, 2), (3, 4)))
const iname_g = joinpath(@__DIR__, "data", "data_20181124_124215.mat")
const params_g, logicals_g = NaCsData.load_striped_mat(iname_g)
data_nacs_g = NaCsData.select_count(params_g, logicals_g, NaCsData.select_single((1, 2), (3, 4)))
const iname_h = joinpath(@__DIR__, "data", "data_20181124_133006.mat")
const params_h, logicals_h = NaCsData.load_striped_mat(iname_h)
data_nacs_h = NaCsData.select_count(params_h, logicals_h, NaCsData.select_single((1, 2), (3, 4)))
const iname_i = joinpath(@__DIR__, "data", "data_20181124_142150.mat")
const params_i, logicals_i = NaCsData.load_striped_mat(iname_i)
data_nacs_i = NaCsData.select_count(params_i, logicals_i, NaCsData.select_single((1, 2), (3, 4)))
const iname_j = joinpath(@__DIR__, "data", "data_20181124_145022.mat")
const params_j, logicals_j = NaCsData.load_striped_mat(iname_j)
data_nacs_j = NaCsData.select_count(params_j, logicals_j, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (298.2055 .+ (-5:0.5:5) .* 1e-3,
                298.2055 .+ (-5:0.5:5) .* 1e-3)
const spec_b = 298.204 .+ (-5:1:5) .* 1e-3
const spec_c = 298.204 .+ (-5:1:5) .* 1e-3
const spec_d = 298.200 .+ (-10:1:10) .* 1e-3
const spec_e = 298.199 .+ (-5:1:5) .* 1e-3
const spec_f = 298.198 .+ (-5:1:5) .* 1e-3
const spec_g = 298.194 .+ (-5:1:5) .* 1e-3
const spec_h = 298.180 .+ (-5:1:5) .* 1e-3
const spec_i = 298.177 .+ (-5:1:5) .* 1e-3
const spec_j = 298.173 .+ (-5:1:5) .* 1e-3

const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)
const split_nacs_b = NaCsData.split_data(data_nacs_b, spec_b)
const split_nacs_c = NaCsData.split_data(data_nacs_c, spec_c)
const split_nacs_d = NaCsData.split_data(data_nacs_d, spec_d)
const split_nacs_e = NaCsData.split_data(data_nacs_e, spec_e)
const split_nacs_f = NaCsData.split_data(data_nacs_f, spec_f)
const split_nacs_g = NaCsData.split_data(data_nacs_g, spec_g)
const split_nacs_h = NaCsData.split_data(data_nacs_h, spec_h)
const split_nacs_i = NaCsData.split_data(data_nacs_i, spec_i)
const split_nacs_j = NaCsData.split_data(data_nacs_j, spec_j)

data_a0 = split_nacs_a[1]
data_a5 = split_nacs_a[2]
data_a10 = split_nacs_b
data_a15 = split_nacs_c
data_a20 = split_nacs_d
data_a25 = split_nacs_e
data_a30 = split_nacs_f
data_a45 = split_nacs_g
data_a60 = split_nacs_h
data_a75 = split_nacs_i
data_a90 = split_nacs_j

const prefix = joinpath(@__DIR__, "imgs", "data_20181124_raman_b_agl")

figure()
NaCsPlot.plot_survival_data(data_a0, fmt="C0o-", label="0\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a5, fmt="C1o-", label="5\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a10, fmt="C2o-", label="10\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a15, fmt="C3o-", label="15\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a20, fmt="C4o-", label="20\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a25, fmt="C5o-", label="25\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a30, fmt="C6o-", label="30\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a45, fmt="C7o-", label="45\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a60, fmt="C8o-", label="60\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a75, fmt="C9o-", label="75\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_a90, fmt="C0s-", label="90\${}^\\circ\$")
grid()
legend(ncol=4, fontsize="small", borderpad=0.3, labelspacing=0.2,
       handletextpad=0.4, columnspacing=0.6, borderaxespad=0.4)
ylim([0, 0.8])
title("Raman spectrum at 60% B field")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_pos")

NaCsPlot.maybe_show()
