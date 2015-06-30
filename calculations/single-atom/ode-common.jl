#!/usr/bin/julia -f

using ODE
using PyPlot

function vector_hcat{T}(V::Vector{Vector{T}})
    height = length(V[1])
    T[V[j][i] for j in 1:length(V), i in 1:height]
end

println("Import done.")
