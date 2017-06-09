#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: interactive
import NaCsCalc: Trap
using Cubature
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 68.8e3, 2π / 589e-9) * 0.67
const η_ax_op = Trap.η(m_Na, 68.8e3, 2π / 589e-9)
const η_ra1 = Trap.η(m_Na, 430e3, 2π / 589e-9)
const η_ra2 = Trap.η(m_Na, 589.5e3, 2π / 589e-9)

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

"""
    op_coupling(n1, n2, η::T, v::T, φ::T, isσ::Bool, cosθ_dri::T, sincosθ_quan::T)

Given an uniform distribution in `[0, 1]` of `v` and an uniform distribution in `[0, π]` of φ
this function will return the distribution of matrix element of matrix element between `n1`
and `n2` motional states after one photon scattering.

The polarization of the emission is determined by `isσ`.
The overlap between the axis and the drive beam direction is `cosθ_dri`.
The `sincos` of the angle between the axis and the quantization axis is `sincosθ_quan`.
"""
function op_coupling(n1, n2, η::T, v::T, φ::T, isσ::Bool, cosθ_dri::T,
                     sincosθ_quan::Tuple{T,T}) where T
    cosθ::T = emission_angle(v, isσ)
    sinθ::T = if abs(cosθ) >= 1
        T(0)
    else
        @fastmath sqrt(1 - cosθ^2)
    end
    sinθ_quan, cosθ_quan = sincosθ_quan
    η_eff = @fastmath η * (cosθ * cosθ_quan + sinθ * sinθ_quan * cos(φ) + cosθ_dri)
    return Trap.sideband(n1, n2, abs(η_eff))^2 / π
end

function op_heating_ax(n1, n2, η::T, isσ) where T
    hcubature(x->op_coupling(n1, n2, η, T(x[1]), T(x[2]), isσ, T(0), (T(1), T(0))),
              [0.0, 0.0], [1.0, π], abstol=1e-6)[1]
end

function op_heating_ra(n1, n2, η::T, isσ) where T
    hcubature(x->op_coupling(n1, n2, η, T(x[1]), T(x[2]), isσ,
                             sqrt(T(0.5)), (sqrt(T(0.5)), sqrt(T(0.5)))),
              [0.0, 0.0], [1.0, π], abstol=1e-6)[1]
end

function op_heating_all(cb::Function, sz1, sz2, η::T, isσ) where T
    res = Matrix{T}(sz1, sz2)
    @inbounds for j in 1:sz2
        for i in 1:sz1
            fill_reverse = false
            if i < j && sz1 > sz2
                continue
            elseif i > j && sz1 <= sz2
                continue
            end
            v = cb(i - 1, j - 1, η, isσ)
            res[i, j] = v
            if i <= sz2 && j <= sz1
                res[j, i] = v
            end
        end
    end
    res
end

const coupling_ax = (op_heating_all(op_heating_ax, 100, 100, Float32(η_ax_op), false) *
                     op_heating_all(op_heating_ax, 100, 100, Float32(η_ax_op), true))
const coupling_ra1 = (op_heating_all(op_heating_ra, 30, 30, Float32(η_ra1), false) *
                      op_heating_all(op_heating_ra, 30, 30, Float32(η_ra1), true))
const coupling_ra2 = (op_heating_all(op_heating_ra, 30, 30, Float32(η_ra2), false) *
                      op_heating_all(op_heating_ra, 30, 30, Float32(η_ra2), true))

function p_heat(coupling, sz)
    T = eltype(coupling)
    res = Vector{T}(sz)
    @inbounds for i in 1:sz
        v = zero(T)
        for j in (i + 1):size(coupling, 1)
            v += coupling[j, i]
        end
        res[i] = v
    end
    return res
end

function maybe_save(name)
    if !interactive()
        dir = dirname(name)
        if !isempty(dir)
            mkpath(dir, 0o755)
        end
        savefig("$name.pdf"; bbox_inches="tight", transparent=true)
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

const nmax = 70
const nmax_r = 15

function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), ".-", label="\$\\Delta n=$(Δn)\$")
    end
end

figure(figsize=[1.5, 1.1] * 4.8)

ax1 = subplot(211)
plot_sidebands(0:nmax, -1:-1:-5, η_ax)
xlim([0, nmax])
ylim([0, 0.75])
grid()
ylabel("\$|\\langle n |e^{ikr}| n + \\Delta n \\rangle|\$")
text(0, 0.65, "\$\\left|\\Delta n\\right|\\!\\!=\\!\\!1\$", color="C0")
text(9, 0.55, "\$\\left|\\Delta n\\right|\\!\\!=\\!\\!2\$", color="C1")
text(21, 0.47, "\$\\left|\\Delta n\\right|\\!\\!=\\!\\!3\$", color="C2")
text(36, 0.44, "\$\\left|\\Delta n\\right|\\!\\!=\\!4\$", color="C3")
text(53, 0.42, "\$\\left|\\Delta n\\right|\\!\\!=\\!5\$", color="C4")
setp(ax1[:get_xticklabels](), visible=false)

ax2 = subplot(212)
subplots_adjust(hspace=0)
plot(0:nmax, p_heat(coupling_ax, nmax + 1), label="Axial")
plot(0:nmax_r, p_heat(coupling_ra1, nmax_r + 1), label="Radial (axis 2)")
plot(0:nmax_r, p_heat(coupling_ra2, nmax_r + 1), label="Radial (axis 3)")
grid()
legend()
ylim([0, 0.57])
xlim([0, nmax])
ylabel("Heating\nprobability")
xlabel("Motional state")

maybe_save(joinpath(@__DIR__, "imgs/raman_op"))

maybe_show()
