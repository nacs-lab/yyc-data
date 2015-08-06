#!/usr/bin/julia -f

import ..Optical: init_phase!, update_phase!

using ..Atomic
using ..Atomic: get_transition_type

@generated function init_phase!{T,NDri}(opt::OpticalCache{T,NDri})
    quote
        trackers = opt.trackers
        $([:(init_phase!(trackers[$i])) for i in 1:NDri]...)
        opt
    end
end

@generated function update_phase!{T,NDri}(opt::OpticalCache{T,NDri}, dt::T)
    quote
        trackers = opt.trackers
        ($([:(update_phase!(trackers[$i], dt)) for i in 1:NDri]...),)
    end
end

@enum SnapshotType SnapshotX SnapshotK

abstract AbstractMeasure

@inline measure_init(::AbstractMeasure, ::Any) = nothing
@inline measure_finalize(::AbstractMeasure, ::Any) = nothing

function measure_snapshot end

immutable DummyMeasure
end

@inline measure_snapshot(::DummyMeasure, args...) = nothing

@enum DecayType DecayNone DecayLeft DecayRight DecayMiddle

# Before the iterations start the measure is called with
#     measure_init(measure, P)
# This should be used to initialize internal states (e.g. buffers).
#
# At each iteration measure is called with
#     measure_snapshot(measure, P, i, ψ, snapshot_type, decay)
# Where
# P is the propagator
# i is the time index 1:(nstep + 1)
# ψ is a snapshot of the current wavefunction
# snapshot_type::SnapshotType is the type of the wavefunction (X or K basis)
# decay::DecayType is the type of the decay (and whether or not it happened)
#
# After the iteration measure_finalize(measure, P) is called to finish up
# the measurement
function propagate{Sys,T}(P::SystemPropagator{Sys,T},
                          ψ0::AbstractMatrix{Complex{T}}, # nele x nstates
                          measure::AbstractMeasure)
    measure_init(measure, P)
    # Disable denormal values
    set_zero_subnormals(true)

    # Manually hoist some load out of the loop to help compiler optimization
    sys = P.sys

    dt = P.dt
    dx = P.dx
    nstep = P.nstep
    nele = P.nele

    motion_cache = P.motion
    E_k = motion_cache.E_k
    P_k = motion_cache.P_k
    E_x = motion_cache.E_x
    P_x2 = motion_cache.P_x2

    optical_cache = P.optical
    coupling_cache = P.coupling

    tmp = P.tmp
    sotmp = P.sotmp

    p_fft! = P.p_fft!
    p_bfft! = P.p_bfft!
    nstates = System.num_states(Sys)

    # Initialization
    init_phase!(optical_cache)

    ψ_norm::T = 0
    @inbounds for j in 1:nstates
        @simd for i in 1:nele
            ψ = ψ0[i, j]
            sotmp[i, j] = ψ
            ψ_norm += abs2(ψ)
        end
    end

    @inbounds for i in 1:(nstep + 1)
        ps_excited = propagate_x1(sys, sotmp, P_x2, 1 / sqrt(ψ_norm), nele)
        if propagate_jump(P, sys, sotmp, nele, ps_excited, dt, measure, i)
            ψ_norm = 1
            continue
        end
    end

    set_zero_subnormals(false)
    measure_finalize(measure, P)
end

@generated function propagate_x1{Sys,N,T}(sys::Sys, sotmp, P_x2::NTuple{N},
                                          ψ_scale::T, nele)
    @meta_expr inline
    P_x2_vars = [gensym(:P_x2) for i in 1:N]
    P_x2_ele_vars = [gensym(:P_x2_ele) for i in 1:N]
    nstates = System.num_states(Sys)
    ψ_vars = [gensym(:ψ) for i in 1:nstates]
    p_vars = [gensym(:p) for i in 1:nstates]

    transitions = System.get_transition_pairs(Sys)
    to_states = Int[]
    for (from, to) in transitions
        to in to_states || push!(to_states, to)
    end

    init_ex = quote
        $([:($(P_x2_vars[i]) = P_x2[$i]) for i in 1:N]...)
        $([:($(p_vars[i])::$T = 0) for i in to_states]...)
    end

    pot_idxs = System.get_potential_idxs(Sys)

    loop_ex = quote
        @simd for j in 1:nele
            # Load the X phase factor. Scale the phase factor since there's
            # less of them compare to the number of states
            $([:($(P_x2_ele_vars[i]) = $(P_x2_vars[i])[j] * ψ_scale)
               for i in 1:N]...)

            # update the wave function
            $([:($(ψ_vars[i]) = sotmp[j, $i] * $(P_x2_ele_vars[pot_idxs[i]]);
                 sotmp[j, $i] = $(ψ_vars[i]))
               for i in 1:nstates]...)

            $([:($(p_vars[i]) += abs2($(ψ_vars[i]))) for i in to_states]...)
        end
    end

    quote
        @meta_expr inline
        @inbounds begin
            $init_ex
            $loop_ex
        end
        ($([p_vars[i] for i in to_states]...),)
    end
end

@generated function propagate_jump{Sys,T}(P, sys::Sys, sotmp, nele, ps_excited,
                                          dt::T, measure, iteration)
    @meta_expr inline
    nstates = System.num_states(Sys)
    transition_pairs = System.get_transition_pairs(Sys)
    ntrans = length(transition_pairs)
    trans_states_idx = Vector{Int}(ntrans)
    to_states = Int[]
    for i in 1:ntrans
        (from, to) = transition_pairs[i]
        to in to_states || push!(to_states, to)
        trans_states_idx[i] = findin(to_states, to)[1]
    end
    p_decay = [gensym(:p) for i in 1:ntrans]

    init_ex = quote
        intern = sys.intern
        transitions = intern.transitions
    end

    # Compute total probability
    prob_ex = quote
        $([:($(p_decay[i]) = (ps_excited[$(trans_states_idx[i])] *
                              transitions[$i].Γ * dt))
           for i in 1:ntrans]...)
        p_decay_total_lin = +($([p_decay[i] for i in 1:ntrans]...))
        p_decay_total = 1 - exp(-p_decay_total_lin)
    end

    # Make dicision
    decide_ex = quote
        p_r = rand(T)
        if p_r >= p_decay_total
            return false
        end
    end

    # Compute thresholds
    decay_ex = quote
        p_r_scale = p_r * p_decay_total_lin / p_decay_total
        p_accum::T = 0
    end

    for i in 1:ntrans
        do_decay = quote
            p_accum += $(p_decay[i])
            if p_accum >= p_r_scale
                p_excited = ps_excited[$(trans_states_idx[i])]
                propagate_do_jump(P, sys, sotmp, nele, p_excited, dt, measure,
                                  iteration, $(Val{i}()))
                return true
            end
        end
        push!(decay_ex.args, do_decay)
    end

    quote
        $init_ex
        $prob_ex
        $decide_ex
        $decay_ex

        # In case there's some rounding error
        return false
    end
end

@generated function propagate_do_jump{Sys,T,TransId}(P, sys::Sys, sotmp, nele,
                                                     p_excited, dt::T, measure,
                                                     iteration,
                                                     _transid::Val{TransId})
    # Giving the arguments a name because of JuliaLang/julia#12474

    Trans = System.get_transition_types(Sys)[TransId]
    # This is a little confusing, the (from, to) of the decay is the
    # opposite of the one for the transition
    To, From = System.get_transition_pairs(Sys)[TransId]

    nstates = System.num_states(Sys)
    ψ_vars = [gensym(:ψ) for i in 1:nstates]

    ax = System.get_quant_axis(Sys)
    cos²θ = (ax * Vec3D(1, 0, 0))^2 / abs2(ax)
    sin²θ = 1 - cos²θ

    init_ex = quote
        decay_dir = rand($T)
        p_fft! = P.p_fft!
        decay_cache = P.optical.decays[$TransId]
        sin_decay = decay_cache.sins
        cos_decay = decay_cache.coss
        tmp = P.tmp
    end

    # Calculate the probabilities
    trans_pol = get_transition_type(Trans)
    # Propability of on axis and off axis emission relative to the
    # quantization axes
    p_ax::T = trans_pol == Trans_π ? 0.1 : 0.2
    p_offax::T = 1 - 2 * p_ax

    # convert to lab axis
    p_ax = p_ax * cos²θ + p_offax / 2 * sin²θ

    # Decide the direction and really do the decay
    do_decay_ex = quote
        ψ_scale::$T = 1 / sqrt(p_excited)
        if decay_dir > $(2 * p_ax)
            # Off axis decay
            decay_type = DecayMiddle
            @inbounds @simd for j in 1:nele
                $([:($(ψ_vars[i]) = Complex{T}(0)) for i in 1:nstates]...)
                $(ψ_vars[To]) = sotmp[j, $From] * ψ_scale
                $([:(sotmp[j, $i] = $(ψ_vars[i])) for i in 1:nstates]...)
            end
        else
            if decay_dir < $p_ax
                ksign = 1
                decay_type = DecayRight
            else
                ksign = -1
                decay_type = DecayLeft
            end
            @inbounds @simd for j in 1:nele
                $([:($(ψ_vars[i]) = Complex{T}(0)) for i in 1:nstates]...)
                $(ψ_vars[To]) = sotmp[j, $From] * ψ_scale
                $(ψ_vars[To]) *= Complex(cos_decay[j], ksign * sin_decay[j])
                $([:(sotmp[j, $i] = $(ψ_vars[i])) for i in 1:nstates]...)
            end
        end
    end

    # Do the post decay measure
    post_measure_ex = quote
        measure_snapshot(measure, P, iteration, sotmp, SnapshotX, decay_type)
        Base.unsafe_copy!(tmp, sotmp)
        p_fft! * tmp
        ψ_scale = 1 / sqrt(T(nele))
        scale!(ψ_scale, tmp)
        measure_snapshot(measure, P, iteration, tmp, SnapshotK, decay_type)
    end

    quote
        $init_ex
        $do_decay_ex
        $post_measure_ex

        nothing
    end
end
