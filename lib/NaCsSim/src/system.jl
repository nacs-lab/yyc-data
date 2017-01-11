#!/usr/bin/julia

module System

import NaCsCalc: Trap
import ..Utils
import ..Samplers
import ..Setup

import ..Utils: WaveFunc

# Atomic state

typealias State{T<:AbstractFloat,N} NTuple{N,WaveFunc{T,3}}

function (::Type{State{T,N}}){T,N}(nx, ny, nz)
    ntuple(x->WaveFunc{T,3}(nx + 1, ny + 1, nz + 1), Val{N})
end

function Utils.zero!(atomic_state::State)
    for ary in atomic_state
        Utils.zero!(ary)
    end
end

function set_ns!{T}(atomic_state::State{T}, ary, nx, ny, nz)
    Utils.zero!(atomic_state)
    push!(ary, (nx + 1, ny + 1, nz + 1), 1)
    return
end

# Initial condition

immutable ThermalInit{Idx,T<:AbstractFloat}
    nx::T
    ny::T
    nz::T
end

function (init::ThermalInit{Idx,T}){Idx,T}(atomic_state::State{T})
    ary = atomic_state[Idx]
    nmax = size(ary)
    nx = Samplers.thermal(init.nx, nmax[1] - 1)
    ny = Samplers.thermal(init.ny, nmax[2] - 1)
    nz = Samplers.thermal(init.nz, nmax[3] - 1)
    set_ns!(atomic_state, ary, nx, ny, nz)
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

function propagate_op!{T,N}(pulse::OPPulse{T}, state::State{T,N}, maxt::T)
    # First, decide which hyperfine state should be pumped and
    # at what time should it happen
    weights = pulse.weights_buf
    for i in 1:N
        weights[i] = abs2(state[i])
    end
    t, n_i = Samplers.decay(pulse.rates, weights)
    0 < t < maxt || return zero(T), true

    # Next, if we want to do a OP, decide which branch it should take
    n_f = Samplers.select(one(T), pulse.branchings[n_i])

    wf = state[n_i]
    # Then, given the initial hyperfine state, pick an initial vibrational state
    v_i = (-).(Samplers.wavefunc(wf), 1)

    # Finally, given the initial hyperfine+vibrational and final hyperfine state,
    # pick the final vibrational state.
    sz_x, sz_y, sz_z = size(wf)
    v_f = Samplers.op(v_i, (sz_x - 1, sz_y - 1, sz_z - 1),
                      pulse.ηs, pulse.ηdri, pulse.isσ[n_f, n_i])
    if v_f[1] < 0
        Utils.zero!(state)
        return zero(T), false
    end
    set_ns!(state, state[n_f], v_f...)
    return maxt - t, true
end

function (pulse::OPPulse{T}){T,N}(state::State{T,N}, extern_state)
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

function (pulse::RamanPulse{T,N1,N2}){T,N1,N2}(state::State{T}, extern_state)
    ary1 = state[N1]
    ary2 = state[N2]
    l1 = length(ary1.sparse)
    l2 = length(ary2.sparse)
    l1 + l2 == 0 && return true
    # For now...
    @assert l1 + l2 == 1
    if l1 == 1
        # forward
        ele = ary1.sparse[1]
        Δn = pulse.Δn
        hf1 = N1
        hf2 = N2
    else
        # backward
        ele = ary2.sparse[1]
        Δn = (-).(pulse.Δn)
        hf1 = N2
        hf2 = N1
    end
    @assert abs2(ele[2]) ≈ 1
    ns_i = ele[1]
    ns_f = (+).(ns_i, Δn)
    if ns_f[1] < 0 || ns_f[2] < 0 || ns_f[3] < 0
        return true
    end
    sz_x, sz_y, sz_z = size(ary1)
    if ns_f[1] >= sz_x || ns_f[2] >= sz_y || ns_f[3] >= sz_z
        Utils.zero!(state)
        return false
    end
    if l1 == 1
        Ω = pulse.Ωs[1][ns_i[1]] * pulse.Ωs[2][ns_i[2]] * pulse.Ωs[3][ns_i[3]]
    else
        Ω = pulse.Ωs[1][ns_f[1]] * pulse.Ωs[2][ns_f[2]] * pulse.Ωs[3][ns_f[3]]
    end
    p = sin(Ω * pulse.t)^2
    if rand() < p
        set_ns!(state, state[hf2], ns_f...)
    end
    return true
end

# External state / measure
immutable HyperFineMeasure
end
function (::HyperFineMeasure){T,N}(state::State{T,N}, extern_state)
    res = Vector{T}(N + 1)
    s = zero(T)
    @inbounds for i in 1:N
        a = abs2(state[i])
        res[i] = a
        s += a
    end
    res[end] = s
    return res
end
Setup.combine_measures(::HyperFineMeasure, m1, m2) = m1 .+ m2
function Setup.finalize_measure(::HyperFineMeasure, m, n)
    len = length(m)
    return (m[1:(len - 1)] ./ m[len], m[len] / n)
end

immutable NBarMeasure
end
function (::NBarMeasure){T,N}(state::State{T,N}, extern_state)
    res = zeros(T, 4)
    for ary in state
        for ele in ary.sparse
            ns, amp = ele
            p = abs2(amp)
            res[1] += (ns[1] - 1) * p
            res[2] += (ns[2] - 1) * p
            res[3] += (ns[3] - 1) * p
            res[4] += p
        end
    end
    return res
end
Setup.combine_measures(::NBarMeasure, m1, m2) = m1 .+ m2
Setup.finalize_measure(::NBarMeasure, m, n) = (m[1:3] ./ m[4], m[4] / n)

end
