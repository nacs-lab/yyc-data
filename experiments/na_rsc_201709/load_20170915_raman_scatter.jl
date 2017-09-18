#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData

# t < 0: F1 down
# t > 0: F2 counter-op
const iname_a = joinpath(@__DIR__, "data", "data_20170915_085614.csv")
# t < 0: F1 up
# t > 0: F2 diagonal
const iname_b = joinpath(@__DIR__, "data", "data_20170915_101807.csv")
# F1 diagonal
const iname_c = joinpath(@__DIR__, "data", "data_20170915_122858.csv")

const data_a = NaCsData.load_count_csv(iname_a)
const data_b = NaCsData.load_count_csv(iname_b)
const data_c = NaCsData.load_count_csv(iname_c)

const f1_down = NaCsData.map_params((i, v)->-v*1e3, data_a[data_a.params .<= 0])
const f2_counterop = NaCsData.map_params((i, v)->v*1e3, data_a[data_a.params .>= 0])

const f1_up = NaCsData.map_params((i, v)->-v*1e3, data_b[data_b.params .<= 0])
const f2_diagonal = NaCsData.map_params((i, v)->v*1e3, data_b[data_b.params .>= 0])

const f1_diagonal = NaCsData.map_params((i, v)->v*1e3, data_c)
