#!/usr/bin/julia -f

import ..Optical: init_phase!, update_phase!

using ..Atomic
using ..Atomic: get_transition_type

@generated function init_phase!{T,NDri}(opt::OpticalCache{T,NDri})
    @meta_expr inline
    quote
        @meta_expr inline
        trackers = opt.trackers
        $([:(init_phase!(trackers[$i])) for i in 1:NDri]...)
        opt
    end
end

@generated function update_phase!{T,NDri}(opt::OpticalCache{T,NDri}, dt::T)
    @meta_expr inline
    quote
        @meta_expr inline
        trackers = opt.trackers
        ($([:(update_phase!(trackers[$i], dt)) for i in 1:NDri]...),)
    end
end

@enum SnapshotType SnapshotX SnapshotK

abstract AbstractMeasure

@inline measure_init(::AbstractMeasure, ::Any) = nothing
@inline measure_finalize(::AbstractMeasure, ::Any) = nothing

function measure_snapshot end

immutable DummyMeasure <: AbstractMeasure
end

@inline measure_snapshot(::DummyMeasure, args...) = nothing

@enum DecayType DecayNone DecayLeft DecayRight DecayMiddle

abstract AbstractSetup

function setup_init end

immutable StaticSetup{Ary} <: AbstractSetup
    ψ0::Ary # nele x nstates
end

@inline function setup_init{Ary,Sys,T}(P::SystemPropagator{Sys,T},
                                       setup::StaticSetup{Ary}, sotmp)
    nstates = System.num_states(Sys)
    nele = P.nele
    ψ_norm::T = 0
    ψ0 = setup.ψ0
    @inbounds for j in 1:nstates
        @simd for i in 1:nele
            ψ = ψ0[i, j]
            sotmp[i, j] = ψ
            ψ_norm += abs2(ψ)
        end
    end
    ψ_norm
end

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
                          setup::AbstractSetup, measure::AbstractMeasure)
    measure_init(measure, P)
    # Disable denormal values
    set_zero_subnormals(true)

    # Manually hoist some load out of the loop to help compiler optimization
    sys = P.sys

    dt = P.dt
    nstep = P.nstep
    nele = P.nele

    motion_cache = P.motion
    P_k = motion_cache.P_k
    P_x2 = motion_cache.P_x2
    P_Es = motion_cache.P_Es
    P_Γs = motion_cache.P_Γs

    optical_cache = P.optical
    coupling_cache = P.coupling

    tmp = P.tmp
    sotmp = P.sotmp

    p_fft! = P.p_fft!
    p_bfft! = P.p_bfft!

    # Initialization
    init_phase!(optical_cache)
    ψ_norm::T = setup_init(P, setup, sotmp)
    inv_sqrt_nele = 1 / sqrt(T(nele))

    @inbounds for i in 1:(nstep + 1)
        update_phase!(optical_cache, dt)
        ps_excited = propagate_x1(sys, sotmp, P_x2, 1 / sqrt(ψ_norm), nele)
        if propagate_jump(P, sys, sotmp, nele, ps_excited, dt, measure, i)
            ψ_norm = 1
            continue
        end
        measure_snapshot(measure, P, i, sotmp, SnapshotX, DecayNone)
        Base.unsafe_copy!(tmp, sotmp)
        p_fft! * tmp
        Base.unsafe_copy!(sotmp, tmp)
        propagate_k(sys, sotmp, P_k, P_Es, inv_sqrt_nele, nele)
        measure_snapshot(measure, P, i, sotmp, SnapshotK, DecayNone)
        Base.unsafe_copy!(tmp, sotmp)
        p_bfft! * tmp
        Base.unsafe_copy!(sotmp, tmp)
        ψ_norm = propagate_x2(sys, sotmp, P_x2, P_Γs, nele)
        propagate_drive(sys, sotmp, nele, optical_cache, coupling_cache)
    end

    set_zero_subnormals(false)
    measure_finalize(measure, P)
end

@generated function propagate_drive{
    Sys,T,N,Idxs}(sys::Sys, sotmp, nele, optical_cache,
                  coupling_cache::CouplingCache{T,N,Idxs})
    @meta_expr inline
    drive_ids = ([drive_id for (drive_id, trans_id) in Idxs]...)
    transition_pairs = System.get_transition_pairs(Sys)
    trans_coupling_pairs = ([transition_pairs[trans_id]
                             for (drive_id, trans_id) in Idxs]...)
    init_ex = quote
        drive_ids = $drive_ids
        trans_coupling_pairs = $trans_coupling_pairs
        trackers = optical_cache.trackers
        drive_trigcaches = optical_cache.drives
        couplings = coupling_cache.couplings
    end
    drive_ex = quote
        @inbounds for k in 1:N
            drive_id = drive_ids[k]
            tracker = trackers[drive_id]
            sins = drive_trigcaches[drive_id].sins
            coss = drive_trigcaches[drive_id].coss
            # Hamiltonian of the spin part is
            # H_σ = Ω (cos(θ_t + θ_x) σ_x + sin(θ_t + θ_x) σ_y)

            # Propagator is
            # P_σ = exp(im H_σ Δt)
            #     = exp(im Ω (cos(θ_t + θ_x) σ_x +
            #                 sin(θ_t + θ_x) σ_y) Δt)
            #     = cos(Ω Δt) + im * (cos(θ_t + θ_x) σ_x +
            #                         sin(θ_t + θ_x) σ_y) * sin(Ω Δt)
            #     = cos(Ω Δt) + im cos(θ_t + θ_x) sin(Ω Δt) σ_x +
            #       im sin(θ_t + θ_x) sin(Ω Δt) σ_y
            #     = [cos(Ω Δt), im exp(im(θ_t + θ_x)) sin(Ω Δt)
            #        im exp(-im(θ_t + θ_x)) sin(Ω Δt), cos(Ω Δt)]

            coupling = couplings[k]
            P_off = coupling.P_off
            P_diag = coupling.P_diag

            # P_σ11 = P_σ22 = cos(Ω Δt)
            # P_σ12 = im exp(im θ_t) exp(im θ_x) sin(Ω Δt)
            # P_σ21 = -P_σ12'

            T11 = P_diag
            T_pre = tracker.exp_t * P_off

            from, to = trans_coupling_pairs[k]

            @simd for j in 1:nele
                # This is a tight loop with complicated operations.
                # It is important to make sure this is as simple as possible
                # so that we don't run out of registers and spill anything
                # to the stack
                ψ_g = sotmp[j, from]
                ψ_e = sotmp[j, to]
                exp_θ_x = Complex(coss[j], sins[j])

                T12 = T_pre * exp_θ_x

                sotmp[j, from] = T11 * ψ_g + T12 * ψ_e
                sotmp[j, to] = T11 * ψ_e - conj(T12) * ψ_g
            end
        end
    end
    quote
        @meta_expr inline
        $init_ex
        $drive_ex

        nothing
    end
end

@generated function propagate_x2{Sys,N,T}(sys::Sys, sotmp,
                                          P_x2::NTuple{N,SoCVector{T}},
                                          P_Γs, nele)
    # We need to apply the phase factors of P_x2 and also the incoherent part
    # of the transformation. Therefore we need to figure out what states can
    # decay and how fast are they decaying.
    @meta_expr inline
    transition_pairs = System.get_transition_pairs(Sys)
    ntrans = length(transition_pairs)
    nstates = System.num_states(Sys)
    to_states = Int[]
    for i in 1:ntrans
        (from, to) = transition_pairs[i]
        to in to_states || push!(to_states, to)
    end
    sort_to_states = sort(to_states)
    none_to_states = Int[]
    for i in 1:nstates
        i in to_states || push!(none_to_states, i)
    end

    init_ex = quote
        ψ_norm::T = 0
    end

    pot_idxs = System.get_potential_idxs(Sys)

    loop_ex = quote
        for i in ($(sort_to_states...),)
            p_x2_single = P_x2[$pot_idxs[i]]
            p_γ = P_Γs[i]
            @simd for j in 1:nele
                ψ = sotmp[j, i] * (p_x2_single[j] * p_γ)
                sotmp[j, i] = ψ
                ψ_norm += abs2(ψ)
            end
        end
        for i in ($(none_to_states...),)
            p_x2_single = P_x2[$pot_idxs[i]]
            @simd for j in 1:nele
                ψ = sotmp[j, i] * p_x2_single[j]
                sotmp[j, i] = ψ
                ψ_norm += abs2(ψ)
            end
        end
    end

    quote
        @meta_expr inline
        @inbounds begin
            $init_ex
            $loop_ex
        end
        ψ_norm
    end
end

@generated function propagate_k{Sys,T}(sys::Sys, sotmp, P_k, P_Es,
                                       ψ_scale::T, nele)
    @meta_expr inline
    nstates = System.num_states(Sys)
    execute_ex = quote
        @inbounds for i in 1:$nstates
            p_E = P_Es[i] * ψ_scale
            @simd for j in 1:nele
                sotmp[j, i] *= P_k[j] * p_E
            end
        end
    end
    quote
        @meta_expr inline
        $execute_ex
        nothing
    end
end

@generated function propagate_x1{Sys,N,T}(sys::Sys, sotmp, P_x2::NTuple{N},
                                          ψ_scale::T, nele)
    @meta_expr inline
    nstates = System.num_states(Sys)
    p_vars = [gensym(:p) for i in 1:nstates]

    transitions = System.get_transition_pairs(Sys)
    to_states = Int[]
    for (from, to) in transitions
        to in to_states || push!(to_states, to)
    end
    sort_to_states = sort(to_states)
    none_to_states = Int[]
    for i in 1:nstates
        i in to_states || push!(none_to_states, i)
    end

    init_ex = quote
        $([:($(p_vars[i])::$T = 0) for i in sort_to_states]...)
    end

    pot_idxs = System.get_potential_idxs(Sys)

    loop_ex = quote
        for i in ($(sort_to_states...),)
            p_x2_single = P_x2[$pot_idxs[i]]
            p_var::$T = 0
            @simd for j in 1:nele
                ψ = sotmp[j, i] * (p_x2_single[j] * ψ_scale)
                sotmp[j, i] = ψ
                p_var += abs2(ψ)
            end
            $([:(i == $i && ($(p_vars[i]) = p_var))
               for i in sort_to_states]...)
        end
        for i in ($(none_to_states...),)
            p_x2_single = P_x2[$pot_idxs[i]]
            @simd for j in 1:nele
                sotmp[j, i] *= p_x2_single[j] * ψ_scale
            end
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
        p_decay_total_lin = +($T(0), $([p_decay[i] for i in 1:ntrans]...))
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
            p_excited = ps_excited[$(trans_states_idx[i])]
            # Filter out small excited state make sure that we don't get
            # accidentally killed by a rounding error or something like that
            if p_excited >= 1e-5 && p_accum >= p_r_scale
                propagate_do_jump(P, sys, sotmp, nele, p_excited, dt, measure,
                                  iteration, $(Val{i}()))
                return true
            end
        end
        push!(decay_ex.args, do_decay)
    end

    quote
        @meta_expr inline
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
                # This unrolled inner loop should be cheap since
                # most of the operations are zeroing
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
                # This unrolled inner loop should be cheap since
                # most of the operations are zeroing
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

abstract MonteCarloMeasure <: AbstractMeasure

function propagate{Sys,T}(P::SystemPropagator{Sys,T}, setup::AbstractSetup,
                          measure::MonteCarloMeasure)
    sub_measure, n = measure_init(measure, P)
    for i in 1:n
        propagate(P, setup, sub_measure)
        measure_snapshot(measure, P, sub_measure)
    end
    measure_finalize(measure, P)
end
