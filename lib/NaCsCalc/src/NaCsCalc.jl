#!/usr/bin/julia

__precompile__(true)

module NaCsCalc

if VERSION >= v"1.0"
    import Statistics: cor, cov, mean, mean!, median, median!, middle,
    quantile, quantile!, std, stdm, var, varm
else
    import Base: cor, cov, mean, mean!, median, median!, middle,
    quantile, quantile!, std, stdm, var, varm
end

export cor, cov, mean, mean!, median, median!, middle,
    quantile, quantile!, std, stdm, var, varm

include("utils.jl")
include("trap.jl")
include("format.jl")
include("atomic.jl")

if VERSION >= v"1.0"
    export linspace
    import .Utils: linspace
end

end
