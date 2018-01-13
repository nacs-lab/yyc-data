#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
import NaCsCalc.Utils: interactive
using NaCsPlot
using PyPlot

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 85.7e3, 2π / 589e-9) * 0.67
const η_rad1 = Trap.η(m_Na, 479e3, 2π / 589e-9) * √(2)
const η_rad2 = Trap.η(m_Na, 492e3, 2π / 589e-9) * √(2)

@show η_ax
@show η_rad1
@show η_rad2

function find_max_mele(η, order)
    mmax = Trap.sideband(0, order, η)
    i = 1
    while true
        m = Trap.sideband(i, i + order, η)
        if m >= mmax
            mmax = m
        else
            break
        end
        i += 1
    end
    return mmax
end

function sideband_max_ratio(η, orders)
    m0 = Trap.sideband(0, 1, η)
    return [m0 / find_max_mele(η, o) for o in orders]
end

@show sideband_max_ratio(η_ax, 1:6)
@show sideband_max_ratio(η_rad1, 1:2)
@show sideband_max_ratio(η_rad2, 1:2)

const prefix = joinpath(@__DIR__, "imgs", "coupling_201801")

const colors = ["k", "r", "y", "g", "c", "b", "m"]
function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), "$(colors[1 - Δn])o-",
             label="\$\\Delta n=$(Δn)\$")
    end
    text(5, 0.8, "\$\\eta\\approx$(@sprintf("%.2f", η))\$", fontsize=30)
    xlim([0, xlim()[2]])
    ylim([0, 1])
    legend(bbox_to_anchor=(1.05, 1), loc=2)
    grid()
    xlabel("n")
    ylabel("\$|\\langle n |e^{ikr}| n + \\Delta n \\rangle|\$")
end

figure()
plot_sidebands(0:90, -6:0, η_ax)
title("Axial Z coupling strength")
NaCsPlot.maybe_save("$(prefix)_z_0.36_0-6")

figure()
plot_sidebands(0:30, -2:0, η_rad1)
title("Radial X coupling strength")
NaCsPlot.maybe_save("$(prefix)_x_0.32_0-2")

figure()
plot_sidebands(0:30, -2:0, η_rad2)
title("Radial Y coupling strength")
NaCsPlot.maybe_save("$(prefix)_y_0.32_0-2")

NaCsPlot.maybe_show()
