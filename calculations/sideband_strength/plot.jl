#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 60e3, 2π / 589e-9) / √(2)
const η_rad = Trap.η(m_Na, 400e3, 2π / 589e-9) * √(2)

const colors = ["k", "r", "y", "g", "c", "b", "m"]
function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), "$(colors[1 - Δn])o-",
             label="\$\\Delta n=$(Δn)\$")
    end
    text(5, 0.8, "\$\\eta\\approx$(@sprintf("%.2f", η))\$", fontsize=30)
    ylim([0, 1])
    legend(bbox_to_anchor=(1.05, 1), loc=2)
    grid()
    xlabel("n")
    ylabel("\$|\\langle n |e^{ikr}| n + \\Delta n \\rangle|\$")
end
figure()
plot_sidebands(0:70, -6:0, η_ax)
title("Axial coupling strength")
savefig(joinpath(ARGS[1], "coupling_0.46_0-6.svg"),
        bbox_inches="tight", transparent=true)
savefig(joinpath(ARGS[1], "coupling_0.46_0-6.png"),
        bbox_inches="tight", transparent=true)
close()
figure()
plot_sidebands(0:25, -2:0, η_rad)
title("Radial coupling strength")
savefig(joinpath(ARGS[1], "coupling_0.35_0-2.svg"),
        bbox_inches="tight", transparent=true)
savefig(joinpath(ARGS[1], "coupling_0.35_0-2.png"),
        bbox_inches="tight", transparent=true)
close()
# show()
