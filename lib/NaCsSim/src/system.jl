#!/usr/bin/julia

module System

import NaCsCalc: Trap
import NaCsCalc.Utils: binomial_estimate
import NaCsCalc.Format: Unc
import ..Samplers
import ..Setup
import ..DecayRabi

# Atomic state

mutable struct StateC
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

struct ThermalInit{Idx,T<:AbstractFloat}
    nx::T
    ny::T
    nz::T
end

function (init::ThermalInit{Idx,T})(state::StateC, rng) where {Idx,T}
    nmax = state.nmax
    nx = Samplers.thermal(init.nx, nmax[1], rng)
    ny = Samplers.thermal(init.ny, nmax[2], rng)
    nz = Samplers.thermal(init.nz, nmax[3], rng)
    set_ns!(state, Idx, nx, ny, nz)
    state.lost = false
    return
end

# Pulses

struct OP{T}
    t::T
    rates::Matrix{T}
    ηs::NTuple{3,T}
    ηdri::NTuple{3,T}
    isσ::Matrix{Bool}
    qax::NTuple{3,T}
    OP{T}(t, rates, ηs, ηdri, isσ, qax=(T(1), T(0), T(0))) where T =
        new(t, rates, ηs, ηdri, isσ, qax)
end

struct OPPulse{T}
    t::T
    rates::Vector{T}
    branchings::Vector{Vector{T}}
    ηs::NTuple{3,T}
    ηdri::NTuple{3,T}
    isσ::Matrix{Bool}
    qtrans::NTuple{3,NTuple{3,T}}
end

OPCache{T} = Tuple{Vector{T},Vector{Vector{T}}}

function compute_cached_op_branching(cache, rates_2d::Matrix{T}) where T
    type_cache = get!(cache, OP{T}) do
        Dict{Matrix{T},OPCache{T}}()
    end::Dict{Matrix{T},OPCache{T}}
    return get!(type_cache, rates_2d) do
        nx, ny = size(rates_2d)
        nx == ny || throw(ArgumentError("Decay rate must be a square matrix"))
        nx >= 1 || throw(ArgumentError("Must have at least one state"))
        rates_1d = Vector{T}(undef, nx)
        branchings = Vector{Vector{T}}(undef, nx)
        @inbounds for i in 1:nx
            r = zero(T)
            @simd for j in 1:nx
                r += rates_2d[j, i]
            end
            rates_1d[i] = r
            if r == 0
                (branchings[i] = zeros(T, nx))[1] = 1
            else
                b = branchings[i] = Vector{T}(undef, nx)
                @simd for j in 1:nx
                    b[j] = rates_2d[j, i] / r
                end
            end
        end
        return rates_1d, branchings
    end
end

function Setup.compile_pulse(pulse::OP{T}, cache) where {T}
    if size(pulse.isσ) != size(pulse.rates)
        throw(ArgumentError("rates and isσ should have the same sizes"))
    end
    rates, branchings = compute_cached_op_branching(cache, pulse.rates)
    return OPPulse{T}(pulse.t, rates, branchings, pulse.ηs, pulse.ηdri, pulse.isσ,
                      Samplers.vec3d_to_trans(pulse.qax))
end

function propagate_op!(pulse::OPPulse{T}, state::StateC, maxt::T, rng) where {T}
    # First, decide which hyperfine state should be pumped and
    # at what time should it happen
    hf0 = state.hf
    v_i = state.n
    nmax = state.nmax

    t = Samplers.decay(pulse.rates[hf0], rng)
    0 < t < maxt || return zero(T), true

    # Next, if we want to do a OP, decide which branch it should take
    hf1 = Samplers.select(one(T), pulse.branchings[hf0], rng)

    # Finally, given the initial hyperfine+vibrational and final hyperfine
    # state, pick the final vibrational state.
    v_f = Samplers.op(v_i, nmax, pulse.ηs, pulse.ηdri, pulse.isσ[hf1, hf0], rng,
                      pulse.qtrans)
    if v_f[1] < 0
        state.lost = true
        return zero(T), false
    end
    set_ns!(state, hf1, v_f...)
    return maxt - t, true
end

function (pulse::OPPulse)(state::StateC, extern_state, rng)
    maxt = pulse.t
    while maxt > 0
        maxt, cont = propagate_op!(pulse, state, maxt, rng)
        cont || return false
    end
    return true
end

struct Raman{T,N1,N2}
    t::T
    Ω::T
    ηs::NTuple{3,T}
    Δn::NTuple{3,Int}
    nmax::NTuple{3,Int}
    Γ::T
    Raman{T,N1,N2}(t, Ω, ηs, Δn, nmax, Γ=0) where {T,N1,N2} = new(t, Ω, ηs, Δn, nmax, Γ)
end

struct RamanPulse{T,N1,N2}
    t::T
    Ω::T
    Δn::NTuple{3,Int}
    Meles::NTuple{3,Vector{T}}
    expΓ::T
end

RamanKey{T} = Tuple{NTuple{3,T},NTuple{3,T},NTuple{3,Int}}
RamanCache{T} = NTuple{3,Vector{T}}

computeMeles(η::T, Δn, nmax) where {T} =
    T[Trap.sideband(n - 1, n - 1 + Δn, η) for n in 1:(nmax + abs(Δn) + 1)]

function compute_cached_raman(cache, ηs::NTuple{3,T}, Δn, nmax) where T
    type_cache = get!(cache, Raman{T}) do
        Dict{RamanKey{T},RamanCache{T}}()
    end::Dict{RamanKey{T},RamanCache{T}}
    return get!(type_cache, (ηs, Δn, nmax)) do
        Meles1 = computeMeles(ηs[1], Δn[1], nmax[1])
        Meles2 = computeMeles(ηs[2], Δn[2], nmax[2])
        Meles3 = computeMeles(ηs[3], Δn[3], nmax[3])
        return Meles1, Meles2, Meles3
    end
end

function Setup.compile_pulse(pulse::Raman{T,N1,N2}, cache) where {T,N1,N2}
    @assert N1 != N2
    Meles = compute_cached_raman(cache, pulse.ηs, pulse.Δn, pulse.nmax)
    return RamanPulse{T,N1,N2}(pulse.t, pulse.Ω, pulse.Δn, Meles, exp(-pulse.Γ * pulse.t))
end

function (pulse::RamanPulse{T,N1,N2})(state::StateC, extern_state, rng) where {T,N1,N2}
    hf0 = state.hf
    v_i = state.n
    nmax = state.nmax

    if hf0 == N1
        # forward
        Δn = pulse.Δn
        hf1 = N2
    elseif hf0 == N2
        # backward
        Δn = .-pulse.Δn
        hf1 = N1
    else
        return true
    end
    v_f = v_i .+ Δn
    if v_f[1] < 0 || v_f[2] < 0 || v_f[3] < 0
        return true
    end
    nmax_x, nmax_y, nmax_z = nmax
    if v_f[1] > nmax_x || v_f[2] > nmax_y || v_f[3] > nmax_z
        state.lost = true
        return false
    end
    if hf0 == N1
        Ω = pulse.Ω * (pulse.Meles[1][v_i[1] + 1] * pulse.Meles[2][v_i[2] + 1] *
                         pulse.Meles[3][v_i[3] + 1])
    else
        Ω = pulse.Ω * (pulse.Meles[1][v_f[1] + 1] * pulse.Meles[2][v_f[2] + 1] *
                         pulse.Meles[3][v_f[3] + 1])
    end
    p = (1 - cos(Ω * pulse.t) * pulse.expΓ) / 2
    if rand(rng) < p
        set_ns!(state, hf1, v_f...)
    end
    return true
end

struct Scatter{T}
    rates::Matrix{T}
    ηs::NTuple{3,T}
    ηdri::NTuple{3,T}
    isσ::Matrix{Bool}
    qax::NTuple{3,T}
    Scatter{T}(rates, ηs, ηdri, isσ, qax=(T(1), T(0), T(0))) where T =
        new(rates, ηs, ηdri, isσ, qax)
end

struct RealRaman{T,N1,N2}
    t::T
    Ω::T
    ηs::NTuple{3,T}
    Δn::NTuple{3,Int}
    nmax::NTuple{3,Int}
    scatters::Vector{Scatter{T}}
    RealRaman{T,N1,N2}(t, Ω, ηs, Δn, nmax, scatters) where {T,N1,N2} =
        new(t, Ω, ηs, Δn, nmax, scatters)
end

struct ScatterPulse{T}
    # From a given state, the probabilities of going into each states
    branchings::Vector{Vector{T}}
    # Max Lamb-Dicke parameter in the three axis for the scattered photon
    ηs::NTuple{3,T}
    # Lamb-Dicke parameter for the drive beam
    ηdri::NTuple{3,T}
    # The polarization of the decay
    isσ::Matrix{Bool}
    qtrans::NTuple{3,NTuple{3,T}}
end

struct RealRamanPulse{T,N1,N2}
    t::T
    Ω::T
    # The change in motional state on different axis
    Δn::NTuple{3,Int}
    # The Raman transition matrix element as a function of the initial state motional level
    Meles::NTuple{3,Vector{T}}
    # Total scattering rates for different HF states
    Γs::Vector{T}
    # Normalized "branching ratio" to each scattering source from different initial HF states
    branchings::Vector{Vector{T}}
    scatters::Vector{ScatterPulse{T}}
end

function num_states(sp::Scatter)
    size1 = size(sp.rates)
    size2 = size(sp.isσ)
    @assert size1 == size2
    @assert size1[1] == size1[2]
    return size1[1]
end

function normalize0!(ary)
    s = sum(ary)
    if s == 0
        return
    end
    ary ./= s
    return
end

function Setup.compile_pulse(pulse::RealRaman{T,N1,N2}, cache) where {T,N1,N2}
    @assert N1 != N2
    Meles = compute_cached_raman(cache, pulse.ηs, pulse.Δn, pulse.nmax)
    ns = length(pulse.scatters)
    if ns == 0
        return RealRamanPulse{T,N1,N2}(pulse.t, pulse.Ω, pulse.Δn, Meles, zeros(max(N1, N2)),
                                       Vector{T}[], ScatterPulse{T}[])
    end
    nhf = num_states(pulse.scatters[1])
    @assert nhf >= N1
    @assert nhf >= N2
    Γs = zeros(T, nhf)
    sc_branchings = [zeros(T, ns) for i in 1:nhf]
    scatters = Vector{ScatterPulse{T}}(undef, ns)
    for i in 1:ns
        st = pulse.scatters[i]
        @assert nhf == num_states(st)
        rates, branchings = compute_cached_op_branching(cache, st.rates)
        for j in 1:nhf
            sc_branchings[j][i] += rates[j]
            Γs[j] += rates[j]
        end
        scatters[i] = ScatterPulse{T}(branchings, st.ηs, st.ηdri, st.isσ,
                                      Samplers.vec3d_to_trans(st.qax))
    end
    for b in sc_branchings
        normalize0!(b)
    end
    return RealRamanPulse{T,N1,N2}(pulse.t, pulse.Ω, pulse.Δn, Meles, Γs, sc_branchings, scatters)
end

function (pulse::RealRamanPulse{T,N1,N2})(state::StateC, extern_state, rng) where {T,N1,N2}
    hf = state.hf
    v = state.n
    nmax = state.nmax
    nmax_x, nmax_y, nmax_z = nmax

    tmax = pulse.t
    while tmax > 0
        if hf == N1
            # forward
            Δn = pulse.Δn
            hf1 = N2
            Γ₁, Γ₂ = pulse.Γs[N1], pulse.Γs[N2]
        elseif hf == N2
            # backward
            Δn = .-pulse.Δn
            hf1 = N1
            Γ₁, Γ₂ = pulse.Γs[N2], pulse.Γs[N1]
        else
            @goto do_op
        end
        v1 = v .+ Δn
        if v1[1] < 0 || v1[2] < 0 || v1[3] < 0
            # Order too high, no Raman
            @goto do_op
        end
        if v1[1] > nmax_x || v1[2] > nmax_y || v1[3] > nmax_z
            # State too high. Lost.
            state.lost = true
            return false
        end
        if hf == N1
            Ω = pulse.Ω * (pulse.Meles[1][v[1] + 1] * pulse.Meles[2][v[2] + 1] *
                             pulse.Meles[3][v[3] + 1])
        else
            Ω = pulse.Ω * (pulse.Meles[1][v1[1] + 1] * pulse.Meles[2][v1[2] + 1] *
                             pulse.Meles[3][v1[3] + 1])
        end
        rabi_param = DecayRabi.Params{T}(Ω, Γ₁, Γ₂)
        tmax2 = tmax
        t, idx, ψ = DecayRabi.propagate_step(rabi_param, tmax, rng)
        tmax -= t
        if idx == 0
            # No decay happened, pick the state and set it
            if rand(rng) < abs2(ψ[1])
                set_ns!(state, hf, v...)
            else
                set_ns!(state, hf1, v1...)
            end
            return true
        end
        if idx == 2
            hf = hf1
            v = v1
        end

        @goto do_scatter
        @label do_op
        if hf > length(pulse.Γs)
            # No scattering on this state, it can stay forever
            break
        end
        t2 = Samplers.decay(pulse.Γs[hf], rng)
        tmax -= t2
        if !(tmax > 0)
            break
        end

        @label do_scatter
        if hf <= length(pulse.branchings)
            # So now we have a scattering even from state `hf` + `v`.
            # First decide which drive it is to blame.
            sidx = Samplers.select(one(T), pulse.branchings[hf], rng)
            scatter = pulse.scatters[sidx]
            # Now figure out the final state
            hf1 = Samplers.select(one(T), scatter.branchings[hf], rng)
            # Finally, figure out the vibrational states
            v = Samplers.op(v, nmax, scatter.ηs, scatter.ηdri, scatter.isσ[hf1, hf], rng,
                            scatter.qtrans)
            hf = hf1
            if v[1] < 0
                state.lost = true
                return false
            end
        end
    end

    set_ns!(state, hf, v...)
    return true
end

struct MultiOP{T}
    t::T
    scatters::Vector{Scatter{T}}
    MultiOP{T}(t, scatters) where {T} = new(t, scatters)
end

struct MultiOPPulse{T}
    t::T
    # Total scattering rates for different HF states
    Γs::Vector{T}
    # Normalized "branching ratio" to each scattering source from different initial HF states
    branchings::Vector{Vector{T}}
    scatters::Vector{ScatterPulse{T}}
end

function Setup.compile_pulse(pulse::MultiOP{T}, cache) where T
    ns = length(pulse.scatters)
    @assert(ns != 0)
    nhf = num_states(pulse.scatters[1])
    Γs = zeros(T, nhf)
    sc_branchings = [zeros(T, ns) for i in 1:nhf]
    scatters = Vector{ScatterPulse{T}}(undef, ns)
    for i in 1:ns
        st = pulse.scatters[i]
        @assert nhf == num_states(st)
        rates, branchings = compute_cached_op_branching(cache, st.rates)
        for j in 1:nhf
            sc_branchings[j][i] += rates[j]
            Γs[j] += rates[j]
        end
        scatters[i] = ScatterPulse{T}(branchings, st.ηs, st.ηdri, st.isσ,
                                      Samplers.vec3d_to_trans(st.qax))
    end
    for b in sc_branchings
        normalize0!(b)
    end
    return MultiOPPulse{T}(pulse.t, Γs, sc_branchings, scatters)
end

function (pulse::MultiOPPulse{T})(state::StateC, extern_state, rng) where {T}
    hf = state.hf
    v = state.n
    nmax = state.nmax

    tmax = pulse.t
    while tmax > 0
        t2 = Samplers.decay(pulse.Γs[hf], rng)
        tmax -= t2
        if !(tmax > 0)
            break
        end

        if hf <= length(pulse.branchings)
            # So now we have a scattering even from state `hf` + `v`.
            # First decide which drive it is to blame.
            sidx = Samplers.select(one(T), pulse.branchings[hf], rng)
            scatter = pulse.scatters[sidx]
            # Now figure out the final state
            hf1 = Samplers.select(one(T), scatter.branchings[hf], rng)
            # Finally, figure out the vibrational states
            v = Samplers.op(v, nmax, scatter.ηs, scatter.ηdri, scatter.isσ[hf1, hf], rng,
                            scatter.qtrans)
            hf = hf1
            if v[1] < 0
                state.lost = true
                return false
            end
        end
    end

    set_ns!(state, hf, v...)
    return true
end

struct Filter{F}
    f::F
end

function (l::Filter)(state::StateC, extern_state, rng)
    hf = state.hf
    v = state.n
    if !l.f(hf, v)
        state.lost = true
        return false
    end
    return true
end

binomial_unc(a, s) = Unc(binomial_estimate(a, s)...)

@inline function check_abort(r, n)
    n < 1000 && return false
    if r * 2 > n
        r = n - r
    end
    p, unc = binomial_estimate(r, n, 1f0)
    return unc < 0.005 || unc < p * 0.01
end

# External state / measure
struct HyperFineMeasure{N}
end
Setup.create_measure(::HyperFineMeasure{N}, seq) where {N} = zeros(Int, N + 1)
function (::HyperFineMeasure{N})(res::Vector{Int}, state::StateC, extern_state, rng) where N
    if !state.lost
        res[state.hf] += 1
        res[N + 1] += 1
    end
    return res
end
function Setup.finalize_measure(::HyperFineMeasure{N}, m, n) where {N}
    total = m[N + 1]
    return (ntuple(i->binomial_unc(m[i], total), Val(N)),
            binomial_unc(total, n))
end
function Setup.abort_measure(::HyperFineMeasure, res::Vector{Int}, n)
    total = res[end]
    total < 1000 && return false
    @inbounds for i in 1:(length(res) - 1)
        check_abort(res[i], total) || return false
    end
    return check_abort(total, n)
end

struct NBarMeasure
end
mutable struct NBarResult
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
function (::NBarMeasure)(res::NBarResult, state::StateC, extern_state, rng)
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
@inline function check_abort_x2(x, x², n)
    xbar² = (x / n)^2
    x²bar = x² / n
    unc² = (x²bar - xbar²) / (n - 1)
    return unc² < 0.000025 / 4 || unc² < xbar² * 0.0001
end
function Setup.abort_measure(::NBarMeasure, res::NBarResult, n)
    n < 1000 && return false
    total = res.n
    check_abort_x2(res.nx, res.nx², total) || return false
    check_abort_x2(res.ny, res.ny², total) || return false
    check_abort_x2(res.nz, res.nz², total) || return false
    return check_abort(total, n)
end

struct FilterMeasure{F}
    cb::F
end
Setup.create_measure(::FilterMeasure, seq) = Ref{Int}(0)
function (measure::FilterMeasure)(res::Ref{Int}, state::StateC, extern_state, rng)
    if !state.lost && measure.cb(state.n, state.hf)
        res[] += 1
    end
    return res
end
Setup.finalize_measure(::FilterMeasure, m, n) = binomial_unc(m[], n)
Setup.abort_measure(::FilterMeasure, m, n) = check_abort(m[], n)

GroundStateMeasure() = FilterMeasure() do n, hf
    n == (0, 0, 0)
end

end
