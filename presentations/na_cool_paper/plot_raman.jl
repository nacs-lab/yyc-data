#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: interactive
import NaCsCalc: Trap
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 67e3, 2π / 589e-9) / √(2)
const η_rad1 = Trap.η(m_Na, 420e3, 2π / 589e-9) * √(2)
const η_rad2 = Trap.η(m_Na, 580e3, 2π / 589e-9) * √(2)

# const colors = ["k", "r", "y", "g", "c", "b", "m"]
const colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), ".-", label="\$\\Delta n=$(Δn)\$")
    end
end

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
plot_sidebands(0:70, -1:-1:-5, η_ax)
xlim([0, 70])
ylim([0, 0.79])
grid()
xlabel("Vibrational state")
ylabel("\$|\\langle n |e^{ikr}| n + \\Delta n \\rangle|\$")
text(0, 0.65, "\$\\left|\\Delta n\\right|\\!\\!=\\!\\!1\$", color="C0")
text(9, 0.55, "\$\\left|\\Delta n\\right|\\!\\!=\\!\\!2\$", color="C1")
text(21, 0.47, "\$\\left|\\Delta n\\right|\\!\\!=\\!\\!3\$", color="C2")
text(36, 0.44, "\$\\left|\\Delta n\\right|\\!\\!=\\!4\$", color="C3")
text(53, 0.42, "\$\\left|\\Delta n\\right|\\!\\!=\\!5\$", color="C4")
maybe_save(joinpath(@__DIR__, "imgs/coupling_ax"))

maybe_show()
