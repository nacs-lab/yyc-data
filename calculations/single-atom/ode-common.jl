#!/usr/bin/julia -f

using ODE
using PyPlot

function vector_hcat{T}(V::Vector{Vector{T}})
    height = length(V[1])
    T[V[j][i] for j in 1:length(V), i in 1:height]
end

# Array only for now
function solve_ode(t0, y0, f, t1, h)
    if t1 <= t0
        error("End time should be after start time.")
    end
    ts = t0:h:t1
    nele = length(y0)
    nstep = length(ts)
    ys = Array{eltype(y0), 2}(nele, nstep)
    ys[:, 1] = y0
    h2 = h / 2
    h3 = h / 3
    h6 = h / 6
    k1 = similar(y0)
    k2 = similar(y0)
    k3 = similar(y0)
    k4 = similar(y0)
    tmp = similar(y0)
    for i in 2:nstep
        t = ts[i]
        prev = sub(ys, (:, i - 1))
        f(t, prev, k1)
        @inbounds for j in 1:nele
            tmp[j] = ys[j, i - 1] + h2 * k1[j]
        end
        f(t + h2, tmp, k2)
        @inbounds for j in 1:nele
            tmp[j] = ys[j, i - 1] + h2 * k2[j]
        end
        f(t + h2, tmp, k3)
        @inbounds for j in 1:nele
            tmp[j] = ys[j, i - 1] + h * k3[j]
        end
        f(t + h, tmp, k4)
        @inbounds for j in 1:nele
            ys[j, i] = ys[j, i - 1] + (h6 * k1[j] + h3 * k2[j] +
                                       h3 * k3[j] + h6 * k4[j])
        end
    end
    collect(ts), ys
end

abstract ODEKernel

function call(k::ODEKernel, t, y)
    ydot = similar(y)
    k(t, y, ydot)
    ydot
end

println("Import done.")
