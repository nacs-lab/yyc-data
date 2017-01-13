#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

const ns = 0:70
const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 60e3, 2π / 589e-9) / √(2)

plot(ns, abs.(Trap.sideband.(ns, ns .- 4, η_ax)), "bo-",
     label="\$\\Delta n=-4\$")
plot(ns, abs.(Trap.sideband.(ns, ns .- 3, η_ax)), "go-",
     label="\$\\Delta n=-3\$")
plot(ns, abs.(Trap.sideband.(ns, ns .- 2, η_ax)), "yo-",
     label="\$\\Delta n=-2\$")
plot(ns, abs.(Trap.sideband.(ns, ns .- 1, η_ax)), "ro-",
     label="\$\\Delta n=-1\$")
text(5, 0.8, "\$\\eta\\approx$(@sprintf("%.2f", η_ax))\$", fontsize=30)
ylim([0, 1])
legend()
grid()
xlabel("n")
ylabel("\$|\\langle n |e^{ikr}| n + \\Delta n \\rangle|\$")
title("Sideband coupling strength")
savefig(joinpath(ARGS[1], "coupling.svg"),
        bbox_inches="tight", transparent=true)
# show()
