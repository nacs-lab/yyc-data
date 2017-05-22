#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
import NaCsCalc.Utils: interactive, sincos
using Cubature
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 67e3, 2π / 589e-9)
const η_ra = Trap.η(m_Na, 580e3, 2π / 589e-9)

function op_coupling_kernel(n1, n2, η::T, θ::T, ϕ::T, isσ) where T
    sinθ, cosθ = sincos(θ)
    @fastmath η_eff = η * cosθ * sin(ϕ)
    factor = if isσ
        T(3 / (2π)) * sinθ^3
    else
        T(3 / (4π)) * (1 + cosθ^2) * sinθ
    end
    return Trap.sideband(n1, n2, η_eff)^2 * factor
end

function op_heating(n1, n2, η::T, isσ) where T
    hcubature(x->op_coupling_kernel(n1, n2, η, T(x[1]), T(x[2]), isσ), [0.0, 0.0], [π / 2, π],
              abstol=1e-5)[1]
end

function op_heating_all(sz1, sz2, η::T, isσ) where T
    res = Matrix{T}(sz1, sz2)
    @inbounds for j in 1:sz2
        for i in 1:sz1
            fill_reverse = false
            if i < j && sz1 > sz2
                continue
            elseif i > j && sz1 <= sz2
                continue
            end
            v = op_heating(i - 1, j - 1, η, isσ)
            res[i, j] = v
            if i <= sz2 && j <= sz1
                res[j, i] = v
            end
        end
    end
    res
end

const coupling_ax = (op_heating_all(100, 100, Float32(η_ax), false) *
                     op_heating_all(100, 100, Float32(η_ax), true))
const coupling_ra = (op_heating_all(100, 100, Float32(η_ra), false) *
                     op_heating_all(100, 100, Float32(η_ra), true))

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

const colors = ["k", "r", "y", "g", "c", "b", "m"]

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

const nmax = 30
const nmax_r = 15
plot(0:nmax, p_heat(coupling_ax, nmax + 1), label="Axial")
plot(0:nmax_r, p_heat(coupling_ra, nmax_r + 1), label="Radial")
grid()
legend()
ylim([0, 0.42])
xlim([0, nmax])
ylabel("Heating probability")
xlabel("Vibrational state")
maybe_save(joinpath(@__DIR__, "imgs/heating_op"))

maybe_show()
