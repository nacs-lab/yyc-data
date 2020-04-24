#!/usr/bin/julia -f

__precompile__(true)

module NaCsData

include("load.jl")
include("selector.jl")
include("fitting.jl")

end
