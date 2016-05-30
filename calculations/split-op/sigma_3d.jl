#!/usr/bin/julia -f

# Simplest unitary system with multiple non-commuting operators: 2 level system

# The system is described by three real numbers, ``c_x``, ``c_y``, ``c_z``
# and the Hamiltonian is ``H = c_x σ_x + c_y σ_y + c_z σ_z``
# We'll evolve the system with exact diagonalization and split operator methods

function exp_n_exact{_T<:Real}(cs::Vector{_T}, dt)
    @assert length(cs) == 3
    T = float(_T)
    CT = Complex{float(T)}
    σ_x = CT[0 1;1 0]
    σ_y = CT[0 -1im;1im 0]
    σ_z = CT[1 0;0 -1]
    σ_n = σ_x * cs[1] + σ_y * cs[2] + σ_z * cs[3]
    return expm(im * T(dt) * σ_n)
end

function exp_n_simple_split{_T<:Real}(cs::Vector{_T}, dt)
    @assert length(cs) == 3
    T = float(_T)
    CT = Complex{T}
    σ_x = CT[0 1;1 0]
    σ_y = CT[0 -1im;1im 0]
    σ_z = CT[1 0;0 -1]
    exp_x = expm(im * T(dt) * cs[1] * σ_x)
    exp_y = expm(im * T(dt) * cs[2] * σ_y)
    exp_z = expm(im * T(dt) * cs[3] * σ_z)
    return exp_x * exp_y * exp_z
end

function propagate{_T<:Real}(cs::Vector{_T}, tmax, dt, _ψ0::Vector, exp_n_getter)
    @assert length(_ψ0) == 2
    T = float(_T)
    CT = Complex{T}
    ψ0 = convert(Vector{CT}, _ψ0)
    nt = round(Int, tmax ÷ dt) + 1
    res = Matrix{CT}(2, nt)
    res[:, 1] = ψ0
    ψ = copy(ψ0)
    exp_n = exp_n_getter(cs, dt)
    for i in 2:nt
        ψ = exp_n * ψ
        res[:, i] = ψ
    end
    return res
end

using PyPlot

function plot_propagate(cs, tmax, dt, ψ0, exp_n_getter, name)
    ts = 0:dt:tmax
    res = propagate(cs, tmax, dt, ψ0, exp_n_getter)
    plot(ts, abs2(res[1, :]), label="$cs $name")
end

plot_exact(cs, tmax, dt, ψ0) =
    plot_propagate(cs, tmax, dt, ψ0, exp_n_exact, "exact")
plot_simple_split(cs, tmax, dt, ψ0) =
    plot_propagate(cs, tmax, dt, ψ0, exp_n_simple_split, "simple")

# plot_exact([1, 1, 0], 10π, 1, [1.0, 0])
# plot_exact([1, 1, 0], 10π, 0.001, [1.0, 0])
plot_exact([1, 1, 0], 10π, 0.1, [1.0, 0])
plot_simple_split([1, 1, 0], 10π, 0.1, [1.0, 0])
# plot_exact([0, 1, 1], 10π, 0.001, [1.0, 0])
# plot_exact([1, 0, 1], 10π, 0.001, [1.0, 0])
legend()
show()
