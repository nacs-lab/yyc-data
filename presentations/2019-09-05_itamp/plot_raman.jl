#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: interactive
import NaCsCalc: Trap
using NaCsPlot
using PyPlot

# matplotlib["rcParams"][:update](Dict("font.weight" => "normal"))

const m_Na = 23e-3 / 6.02e23
const η_rx = Trap.η(m_Na, 479e3, 2π / 589e-9) * √(2) * 0.96
const η_ry = Trap.η(m_Na, 492e3, 2π / 589e-9) * √(2) * 0.933
const η_az = Trap.η(m_Na, 85.7e3, 2π / 589e-9) * 0.67
const η_az_op = Trap.η(m_Na, 85.7e3, 2π / 589e-9)

"""
    emission(v, isσ::Bool) -> (sinθ, cosθ)

Given an uniform distribution in `[0, 1]` of `v` this function will return a distribution
of `cosθ` that matches a dipole emission pattern (the type of the transition is determined by
`isσ`).
"""
function emission_angle(v::T, isσ::Bool) where T<:AbstractFloat
    # Returns `cosθ`. The caller should be ready to handle `|cosθ| > 1`
    # The PDFs of the `θ` distribution are
    # `3 / 4 * (1 - cos²θ) * sinθ` for π light and
    # `3 / 8 * (1 + cos²θ) * sinθ` for σ± light
    # The corresponding CDFs are
    # `1 / 4 * (3cosθ - cos³θ) + 1 / 2` for π light and
    # `1 / 8 * (3cosθ + cos³θ) + 1 / 2` for σ± light
    # For a random number `v` we need to solve
    # `3cosθ - cos³θ = 4v - 2` for π light and
    # `3cosθ + cos³θ = 8v - 4` for σ± light
    # The (real) solution is
    # `cosθ = x + 1 / x` where `x = ∛(2√(v² - v) - 2v + 1)` for π light and
    # `cosθ = x - 1 / x` where `x = ∛(√(16v² - 16v + 5) + 4v - 2)` for σ± light
    # Note that the `x` for π polarization is complex
    if isσ
        y = muladd(v, 4, -2)
        x = @fastmath cbrt(sqrt(muladd(y, y, 1)) + y)
        return x - 1 / x
    else
        θ′ = @fastmath acos(muladd(T(-2), v, T(1)))
        cosθ = @fastmath cos(muladd(θ′, T(1 / 3), - T(2π / 3))) * 2
        return cosθ
    end
end

const nmax = 95
const nmax_r = 15

function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), ".-", label="\$\\Delta n=$(Δn)\$")
    end
end

figure()
plot_sidebands(0:nmax, -1:-1, η_az)
xlim([0, nmax])
ylim([0, 0.75])
grid()
ylabel("\$|\\langle n |e^{i\\Delta\\vec k\\cdot\\vec z}| n + \\Delta n \\rangle|\$")
text(0, 0.61, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!1\$", color="C0")
xlabel("Motional state \$n\$")
NaCsPlot.maybe_save(joinpath(@__DIR__, "raman_rabi_1"))

figure()
plot_sidebands(0:nmax, -1:-1:-5, η_az)
xlim([0, nmax])
ylim([0, 0.75])
grid()
ylabel("\$|\\langle n |e^{i\\Delta\\vec k\\cdot\\vec z}| n + \\Delta n \\rangle|\$")
text(0, 0.61, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!1\$", color="C0")
text(12, 0.51, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!2\$", color="C1")
text(29, 0.46, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!3\$", color="C2")
text(47, 0.42, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!4\$", color="C3")
text(74, 0.40, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!5\$", color="C4")
xlabel("Motional state \$n\$")
NaCsPlot.maybe_save(joinpath(@__DIR__, "raman_rabi_all"))

NaCsPlot.maybe_show()
