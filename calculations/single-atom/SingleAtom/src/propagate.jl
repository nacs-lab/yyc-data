#!/usr/bin/julia -f

import ..Optical: init_phase!, update_phase!

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
        propagate_x1(sys, sotmp, P_x2, 1 / sqrt(ψ_norm), nele)
    end

    set_zero_subnormals(false)
    measure_finalize(measure, P)
end

@generated function propagate_x1{Sys,N,T}(sys::Sys, sotmp,
                                          P_x2::NTuple{N}, ψ_scale::T, nele)
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
