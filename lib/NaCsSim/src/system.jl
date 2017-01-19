#!/usr/bin/julia

module System

import NaCsCalc: Trap
import NaCsCalc.Utils: binomial_estimate
import NaCsCalc.Format: Unc
import ..Samplers
import ..Setup

# Atomic state

type StateC
    nmax::NTuple{3,Int}
    n::NTuple{3,Int}
    hf::Int
    lost::Bool
    function StateC(nx, ny, nz)
        new((nx, ny, nz), (0, 0, 0), 1, false)
    end
end

@inline function set_ns!(state::StateC, hf, nx, ny, nz)
    state.n = (nx, ny, nz)
    state.hf = hf
    return
end

# Initial condition

immutable ThermalInit{Idx,T<:AbstractFloat}
    nx::T
    ny::T
    nz::T
end

function (init::ThermalInit{Idx,T}){Idx,T}(state::StateC)
    nmax = state.nmax
    nx = Samplers.thermal(init.nx, nmax[1])
    ny = Samplers.thermal(init.ny, nmax[2])
    nz = Samplers.thermal(init.nz, nmax[3])
    set_ns!(state, Idx, nx, ny, nz)
    state.lost = false
    return
end

# Pulses

immutable OP{T}
    t::T
    rates::Matrix{T}
    ηs::NTuple{3,T}
    ηdri::NTuple{3,T}
    isσ::Matrix{Bool}
end

immutable OPPulse{T}
    t::T
    rates::Vector{T}
    branchings::Vector{Vector{T}}
    weights_buf::Vector{T}
    ηs::NTuple{3,T}
    ηdri::NTuple{3,T}
    isσ::Matrix{Bool}
end

typealias OPCache{T} Tuple{Vector{T},Vector{Vector{T}},Vector{T}}

function Setup.compile_pulse{T}(pulse::OP{T}, cache)
    if size(pulse.isσ) != size(pulse.rates)
        throw(ArgumentError("rates and isσ should have the same sizes"))
    end
    type_cache = get!(cache, OP{T}) do
        Dict{Matrix{T},OPCache{T}}()
    end::Dict{Matrix{T},OPCache{T}}
    rates, branchings, weights_buf = get!(type_cache, pulse.rates) do
        local rates_2d, rates_1d, branchings, weights_buf
        rates_2d = pulse.rates
        nx, ny = size(rates_2d)
        nx == ny || throw(ArgumentError("Decay rate must be a square matrix"))
        nx >= 1 || throw(ArgumentError("Must have at least one state"))
        rates_1d = Vector{T}(nx)
        branchings = Vector{Vector{T}}(nx)
        weights_buf = Vector{T}(nx)
        @inbounds for i in 1:nx
            r = zero(T)
            @simd for j in 1:nx
                r += rates_2d[j, i]
            end
            rates_1d[i] = r
            if r == 0
                (branchings[i] = zeros(T, nx))[1] = 1
            else
                b = branchings[i] = Vector{T}(nx)
                @simd for j in 1:nx
                    b[j] = rates_2d[j, i] / r
                end
            end
        end
        return rates_1d, branchings, weights_buf
    end
    return OPPulse{T}(pulse.t, rates, branchings, weights_buf,
                      pulse.ηs, pulse.ηdri, pulse.isσ)
end

function propagate_op!{T}(pulse::OPPulse{T}, state::StateC, maxt::T)
    # First, decide which hyperfine state should be pumped and
    # at what time should it happen
    hf0 = state.hf
    v_i = state.n
    nmax = state.nmax

    t = Samplers.decay(pulse.rates[hf0])
    0 < t < maxt || return zero(T), true

    # Next, if we want to do a OP, decide which branch it should take
    hf1 = Samplers.select(one(T), pulse.branchings[hf0])

    # Finally, given the initial hyperfine+vibrational and final hyperfine
    # state, pick the final vibrational state.
    v_f = Samplers.op(v_i, nmax, pulse.ηs, pulse.ηdri, pulse.isσ[hf1, hf0])
    if v_f[1] < 0
        state.lost = true
        return zero(T), false
    end
    set_ns!(state, hf1, v_f...)
    return maxt - t, true
end

function (pulse::OPPulse)(state::StateC, extern_state)
    maxt = pulse.t
    while maxt > 0
        maxt, cont = propagate_op!(pulse, state, maxt)
        cont || return false
    end
    return true
end

immutable Raman{T,N1,N2}
    t::T
    Ω::T
    ηs::NTuple{3,T}
    Δn::NTuple{3,Int}
    nmax::NTuple{3,Int}
end

immutable RamanPulse{T,N1,N2}
    t::T
    Δn::NTuple{3,Int}
    Ωs::NTuple{3,Vector{T}}
end

typealias RamanKey{T} Tuple{NTuple{3,T},NTuple{3,T},NTuple{3,Int}}
typealias RamanCache{T} NTuple{3,Vector{T}}

computeΩs{T}(Ω::T, η::T, Δn, nmax) =
    T[Trap.sideband(n - 1, n - 1 + Δn, η) * Ω for n in 1:(nmax + abs(Δn) + 1)]

function Setup.compile_pulse{T,N1,N2}(pulse::Raman{T,N1,N2}, cache)
    @assert N1 != N2
    type_cache = get!(cache, Raman{T}) do
        Dict{RamanKey{T},RamanCache{T}}()
    end::Dict{RamanKey{T},RamanCache{T}}
    Ωs = get!(type_cache, (pulse.ηs, pulse.Δn, pulse.nmax)) do
        Ωs1 = computeΩs(pulse.Ω, pulse.ηs[1], pulse.Δn[1], pulse.nmax[1])
        Ωs2 = computeΩs(pulse.Ω, pulse.ηs[2], pulse.Δn[2], pulse.nmax[2])
        Ωs3 = computeΩs(pulse.Ω, pulse.ηs[3], pulse.Δn[3], pulse.nmax[3])
        return Ωs1, Ωs2, Ωs3
    end
    return RamanPulse{T,N1,N2}(pulse.t, pulse.Δn, Ωs)
end

function (pulse::RamanPulse{T,N1,N2}){T,N1,N2}(state::StateC, extern_state)
    hf0 = state.hf
    v_i = state.n
    nmax = state.nmax

    if hf0 == N1
        # forward
        Δn = pulse.Δn
        hf1 = N2
    elseif hf0 == N2
        # backward
        Δn = (-).(pulse.Δn)
        hf1 = N1
    else
        return true
    end
    v_f = (+).(v_i, Δn)
    if v_f[1] < 0 || v_f[2] < 0 || v_f[3] < 0
        return true
    end
    nmax_x, nmax_y, nmax_z = nmax
    if v_f[1] > nmax_x || v_f[2] > nmax_y || v_f[3] > nmax_z
        state.lost = true
        return false
    end
    if hf0 == N1
        Ω = (pulse.Ωs[1][v_i[1] + 1] * pulse.Ωs[2][v_i[2] + 1] *
              pulse.Ωs[3][v_i[3] + 1])
    else
        Ω = (pulse.Ωs[1][v_f[1] + 1] * pulse.Ωs[2][v_f[2] + 1] *
              pulse.Ωs[3][v_f[3] + 1])
    end
    p = sin(Ω * pulse.t)^2
    if rand() < p
        set_ns!(state, hf1, v_f...)
    end
    return true
end

binomial_unc(a, s) = Unc(binomial_estimate(a, s)...)

# External state / measure
immutable HyperFineMeasure{N}
end
Setup.create_measure{N}(::HyperFineMeasure{N}, seq) = zeros(Int, N + 1)
function (::HyperFineMeasure{N}){N}(res::Vector{Int}, state::StateC,
                                    extern_state)
    if !state.lost
        res[state.hf] += 1
        res[N + 1] += 1
    end
    return res
end
function Setup.finalize_measure{N}(::HyperFineMeasure{N}, m, n)
    total = m[N + 1]
    return (ntuple(i->binomial_unc(m[i], total), Val{N}),
            binomial_unc(total, n))
end

immutable NBarMeasure
end
type NBarResult
    nx::Int
    ny::Int
    nz::Int
    n::Int
    nx²::Float64
    ny²::Float64
    nz²::Float64
    NBarResult() = new(0, 0, 0, 0,
                       0.0, 0.0, 0.0)
end
Setup.create_measure(::NBarMeasure, seq) = NBarResult()
function (::NBarMeasure)(res::NBarResult, state::StateC, extern_state)
    state.lost && return res
    n = state.n
    res.nx += n[1]
    res.ny += n[2]
    res.nz += n[3]
    res.n += 1
    res.nx² += n[1]^2
    res.ny² += n[2]^2
    res.nz² += n[3]^2
    return res
end
function Setup.finalize_measure(::NBarMeasure, res::NBarResult, n)
    total = res.n
    nx = res.nx / total
    ny = res.ny / total
    nz = res.nz / total
    nx² = res.nx² / total
    ny² = res.ny² / total
    nz² = res.nz² / total
    factor = 1 / sqrt(total - 1)
    σnx = sqrt(nx² - nx^2) * factor
    σny = sqrt(ny² - ny^2) * factor
    σnz = sqrt(nz² - nz^2) * factor
    ((Unc(nx, σnx), Unc(ny, σny), Unc(nz, σnz)), binomial_unc(total, n))
end

immutable FilterMeasure{F}
    cb::F
end
Setup.create_measure(::FilterMeasure, seq) = Ref{Int}(0)
function (measure::FilterMeasure)(res::Ref{Int}, state::StateC, extern_state)
    if !state.lost && measure.cb(state.n, state.hf)
        res[] += 1
    end
    return res
end
Setup.finalize_measure(::FilterMeasure, m, n) = binomial_unc(m[], n)

GroundStateMeasure() = FilterMeasure() do n, hf
    n == (0, 0, 0)
end

end
