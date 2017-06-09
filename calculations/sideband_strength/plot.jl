#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
import NaCsCalc.Utils: interactive
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 68.8e3, 2π / 589e-9) * 0.67
const η_rad1 = Trap.η(m_Na, 430e3, 2π / 589e-9) * √(2)
const η_rad2 = Trap.η(m_Na, 589.5e3, 2π / 589e-9) * √(2)

function calc_fraction(η)
    m1 = Trap.sideband(0, 1, η)
    m2 = Trap.sideband(1, 2, η)
    @assert m1 < m2
    mmax = m2
    i = 2
    while true
        m = Trap.sideband(i, i + 1, η)
        if m > mmax
            mmax = m
        else
            break
        end
    end
    v1 = sin(m1 / mmax * π / 2)^2
    v2 = sin(m2 / mmax * π / 2)^2
    return v2 / v1
end

@show calc_fraction(η_ax)
@show calc_fraction(η_rad1)
@show calc_fraction(η_rad2)

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

const save_fig = get(ENV, "NACS_SAVE_FIG", "true") == "true"

function maybe_save(name)
    if !interactive()
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

figure()
plot_sidebands(0:70, -6:0, η_ax)
title("Axial coupling strength")
maybe_save(joinpath(@__DIR__, "imgs/coupling_0.43_0-6"))

figure()
plot_sidebands(0:25, -2:0, η_rad1)
title("Radial 2 coupling strength")
maybe_save(joinpath(@__DIR__, "imgs/coupling_0.35_0-2"))

figure()
plot_sidebands(0:25, -2:0, η_rad2)
title("Radial 3 coupling strength")
maybe_save(joinpath(@__DIR__, "imgs/coupling_0.29_0-2"))

maybe_show()
