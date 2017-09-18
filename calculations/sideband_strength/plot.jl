#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
import NaCsCalc.Utils: interactive
using NaCsPlot
using PyPlot

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 87.5e3, 2π / 589e-9) * 0.67
const η_rad1 = Trap.η(m_Na, 478e3, 2π / 589e-9) * √(2)
const η_rad2 = Trap.η(m_Na, 499e3, 2π / 589e-9) * √(2)

# function calc_fraction(η)
#     m1 = Trap.sideband(0, 1, η)
#     m2 = Trap.sideband(1, 2, η)
#     @assert m1 < m2
#     mmax = m2
#     i = 2
#     while true
#         m = Trap.sideband(i, i + 1, η)
#         if m > mmax
#             mmax = m
#         else
#             break
#         end
#         i += 1
#     end
#     v1 = sin(m1 / mmax * π / 2)^2
#     v2 = sin(m2 / mmax * π / 2)^2
#     return v2 / v1
# end

# @show calc_fraction(η_ax)
# @show calc_fraction(η_rad1)
# @show calc_fraction(η_rad2)

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
plot_sidebands(0:70, -6:0, η_ax)
title("Axial Z coupling strength")
NaCsPlot.maybe_save(joinpath(@__DIR__, "imgs/coupling_z_0.36_0-6"))

figure()
plot_sidebands(0:25, -2:0, η_rad1)
title("Radial X coupling strength")
NaCsPlot.maybe_save(joinpath(@__DIR__, "imgs/coupling_x_0.32_0-2"))

figure()
plot_sidebands(0:25, -2:0, η_rad2)
title("Radial Y coupling strength")
NaCsPlot.maybe_save(joinpath(@__DIR__, "imgs/coupling_y_0.32_0-2"))

NaCsPlot.maybe_show()
