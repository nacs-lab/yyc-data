#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
using Cubature
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 67e3, 2π / 589e-9)

function op_coupling_kernel(n1, n2, η, θ, ϕ, isσ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    η_eff = η * cosθ * sin(ϕ)
    factor = if isσ
        (3 / (2π)) * sinθ^3
    else
        (3 / (4π)) * (1 + cosθ^2) * sinθ
    end
    return Trap.sideband(n1, n2, η_eff)^2 * factor
end

function op_heating(n1, n2, η, isσ)
    hcubature(x->op_coupling_kernel(n1, n2, η, x[1], x[2], isσ), [0.0, 0.0], [π / 2, π],
              abstol=1e-5)[1]
end

const colors = ["k", "r", "y", "g", "c", "b", "m"]
const save_fig = get(ENV, "NACS_SAVE_FIG", "true") == "true"

function maybe_save(name)
    if save_fig
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if !save_fig
        show()
    end
end

function plot_sidebands(n1s, n2s, η, isσ)
    for i in 1:length(n2s)
        n2 = n2s[i]
        plot(n1s, op_heating.(n1s, n2, η, isσ), "$(colors[i])o-", label="\$n_{init}=$n2\$")
    end
    text(10, 0.8, "\$\\eta\\approx$(@sprintf("%.2f", η))\$", fontsize=30)
    ylim([0, 1])
    xlim([0, xlim()[2]])
    legend(bbox_to_anchor=(1.05, 1), loc=2)
    grid()
    xlabel("\$n_{final}\$")
    ylabel("Probability")
end


figure()
plot_sidebands(0:60, [0, 1, 2, 5, 10, 20, 50], η_ax, false)
title("Axial OP vibrational state distribution\nafter \$\\pi\$ emission")
maybe_save(joinpath(@__DIR__, "imgs/coupling_0.61_op-pi"))

figure()
plot_sidebands(0:60, [0, 1, 2, 5, 10, 20, 50], η_ax, true)
title("Axial OP vibrational state distribution\nafter \$\\sigma\$ emission")
maybe_save(joinpath(@__DIR__, "imgs/coupling_0.61_op-sigma"))

maybe_show()
