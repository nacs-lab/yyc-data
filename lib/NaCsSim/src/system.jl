#!/usr/bin/julia

module System

import ..Utils
import ..Samplers
import ..Setup

import ..Utils: HybridArray

# Atomic state

typealias State{T<:AbstractFloat,N} NTuple{N,HybridArray{T,3}}

function (::Type{State{T,N}}){T,N}(nx, ny, nz)
    ntuple(x->HybridArray{T,3}(nx + 1, ny + 1, nz + 1), Val{N})
end

function Utils.zero!(atomic_state::State)
    for ary in atomic_state
        Utils.zero!(ary)
    end
end

function set_ns!{T}(atomic_state::State{T}, ary, nx, ny, nz)
    Utils.zero!(atomic_state)
    ary.sum = 1
    push!(ary.sparse, ((nx + 1, ny + 1, nz + 1), Complex{T}(1)))
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
    nmax = size(ary.full)
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
        weights[i] = state[i].sum
    end
    t, n_i = Samplers.decay(pulse.rates, weights)
    0 < t < maxt || return zero(T), true

    # Next, if we want to do a OP, decide which branch it should take
    n_f = Samplers.select(one(T), pulse.branchings[n_i])

    wf = state[n_i]
    # Then, given the initial hyperfine state, pick an initial vibrational state
    v_i = Samplers.wavefunction(wf)

    # Finally, given the initial hyperfine+vibrational and final hyperfine state,
    # pick the final vibrational state.
    sz_x, sz_y, sz_z = size(wf.full)
    v_f = Samplers.op(v_i, (sz_x - 1, sz_y - 1, sz_z - 1),
                      pulse.ηs, pulse.ηdri, pulse.isσ[n_f, n_i])
    if v_f[1] < 0
        Utils.zero!(state)
        return zero(T), false
    end
    set_ns!(state, wf, v_f...)
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

# External state / measure
immutable HyperFineMeasure
end
(::HyperFineMeasure)(state::State, extern_state) = [a.sum for a in state]

end
