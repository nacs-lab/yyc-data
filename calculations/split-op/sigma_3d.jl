#!/usr/bin/julia -f

# Simplest unitary system with multiple non-commuting operators: 2 level system

# The system is described by three real numbers, ``c_x``, ``c_y``, ``c_z``
# and the Hamiltonian is ``H = c_x σ_x + c_y σ_y + c_z σ_z``
# We'll evolve the system with exact diagonalization and split operator methods

function get_σ{T<:Real}(cs::Vector{T})
    @assert length(cs) == 3
    CT = Complex{float(T)}
    σ_x = CT[0 1;1 0]
    σ_y = CT[0 -1im;1im 0]
    σ_z = CT[1 0;0 -1]
    return σ_x * cs[1] + σ_y * cs[2] + σ_z * cs[3]
end

function propagate_exact{_T<:Real}(cs::Vector{_T}, tmax, dt, _ψ0::Vector)
    ts = 0:dt:tmax
    @assert length(_ψ0) == 2
    T = float(_T)
    CT = Complex{T}
    ψ0 = convert(Vector{CT}, _ψ0)
    σ_n = get_σ(cs)
    nt = length(ts)
    res = Matrix{CT}(2, nt)
    for i in 1:nt
        res[:, i] = expm(im * T(ts[i]) * σ_n) * ψ0
    end
    return res
end

function propagate_simple_split{_T<:Real}(cs::Vector{_T}, tmax, dt, _ψ0::Vector)
    @assert length(_ψ0) == 2
    T = float(_T)
    CT = Complex{T}
    ψ0 = convert(Vector{CT}, _ψ0)
    nt = round(Int, tmax ÷ dt) + 1
    σ_x = CT[0 1;1 0]
    σ_y = CT[0 -1im;1im 0]
    σ_z = CT[1 0;0 -1]
    res = Matrix{CT}(2, nt)
    res[:, 1] = ψ0
    ψ = copy(ψ0)

    exp_x = expm(im * T(dt) * cs[1] * σ_x)
    exp_y = expm(im * T(dt) * cs[2] * σ_y)
    exp_z = expm(im * T(dt) * cs[3] * σ_z)
    exp_n = exp_x * exp_y * exp_z
    # For precision only so using the simple and slow approach...
    for i in 2:nt
        ψ = exp_n * ψ
        res[:, i] = ψ
    end
    return res
end

using PyPlot

function plot_exact(cs, tmax, dt, ψ0)
    ts = 0:dt:tmax
    res = propagate_exact(cs, tmax, dt, ψ0)
    plot(ts, abs2(res[1, :]), label="$cs exact")
end

function plot_simple_split(cs, tmax, dt, ψ0)
    ts = 0:dt:tmax
    res = propagate_simple_split(cs, tmax, dt, ψ0)
    plot(ts, abs2(res[1, :]), label="$cs simple")
end

plot_exact([1, 1, 0], 10π, 0.1, [1.0, 0])
plot_simple_split([1, 1, 0], 10π, 0.1, [1.0, 0])
# plot_exact([0, 1, 1], 10π, 0.001, [1.0, 0])
# plot_exact([1, 0, 1], 10π, 0.001, [1.0, 0])
legend()
show()
