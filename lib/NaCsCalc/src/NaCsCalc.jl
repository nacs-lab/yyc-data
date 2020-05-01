#!/usr/bin/julia

__precompile__(true)

module NaCsCalc

import Statistics: cor, cov, mean, mean!, median, median!, middle,
    quantile, quantile!, std, stdm, var, varm
export cor, cov, mean, mean!, median, median!, middle,
    quantile, quantile!, std, stdm, var, varm

include("utils.jl")
include("trap.jl")
include("format.jl")
include("atomic.jl")

export linspace
import .Utils: linspace

end
