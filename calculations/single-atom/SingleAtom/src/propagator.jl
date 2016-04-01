#!/usr/bin/julia -f

module Propagate

using Compat
using ..Utils
using ..System
using ..Atomic
import ..Optical

export SystemPropagator, AbstractMeasure, DummyMeasure, MonteCarloMeasure
export AbstractSetup, StaticSetup, propagate

export SnapshotType, SnapshotX, SnapshotK
export DecayType, DecayNone, DecayLeft, DecayRight, DecayMiddle

immutable HMotionCache{T,N,NState}
    E_k::Vector{T}
    P_k::SoCVector{T}

    E_x::NTuple{N,Vector{T}}
    P_x2::NTuple{N,SoCVector{T}}

    P_Es::NTuple{NState,Complex{T}}
    P_Γs::NTuple{NState,T}
end

@generated function HMotionCache{M<:MotionSystem}(sys::M, _dx, _dt, nele)
    nstates = System.num_states(M)
    NPots = length(System.get_potential_types(M))
    T = System.get_value_type(M)
    init_expr = quote
        dx = $T(_dx)
        dt = $T(_dt)
        m = sys.mass
        k0 = $T(2π) / (nele * dx)
        nele_2 = nele ÷ 2
        E_k = Vector{$T}(nele)
        P_k = StructOfArrays(Complex{$T}, nele)
        E_x = ($([:(Vector{$T}(nele)) for i in 1:NPots]...),)
        P_x2 = ($([:(StructOfArrays(Complex{$T}, nele)) for i in 1:NPots]...),)
    end
    calc_k = quote
        @inbounds for i in 1:nele
            k1 = i - 1
            k2 = k1 - nele
            k = ifelse(k1 + k2 <= 0, k1, k2) * k0

            e_k = get_kinetic(m, k)
            E_k[i] = e_k
            P_k[i] = exp((1im) * e_k * dt)
        end
    end
    e_temps = [gensym() for i in 1:NPots]
    calc_x = quote
        @inbounds for i in 1:nele
            x = (i - nele_2) * dx
            $([:($(e_temps[i]) = get_potential(sys.potentials[$i], m, x))
               for i in 1:NPots]...)
            $([:(E_x[$i][i] = $(e_temps[i])) for i in 1:NPots]...)
            $([:(P_x2[$i][i] = exp((1im) * dt / 2 * $(e_temps[i])))
               for i in 1:NPots]...)
        end
    end
    calc_es = quote
        energies = sys.intern.energies
        P_Es = ($([:(exp((1im) * dt * energies[$i])) for i in 1:nstates]...),)
    end
    transition_pairs = System.get_transition_pairs(M)
    ntrans = length(transition_pairs)
    nstates = System.num_states(M)
    decay_rate = Any[[:($T(0))] for i in 1:nstates]
    for i in 1:ntrans
        (from, to) = transition_pairs[i]
        push!(decay_rate[to], :(transitions[$i].Γ))
    end
    calc_γs = quote
        transitions = sys.intern.transitions
        P_Γs = ($([:(exp(-dt / 2 * +($(decay_rate[i]...))))
                    for i in 1:nstates]...),)
    end
    quote
        $init_expr
        $calc_k
        $calc_x
        $calc_es
        $calc_γs
        HMotionCache{$T,$NPots,$nstates}(E_k, P_k, E_x, P_x2, P_Es, P_Γs)
    end
end

# Optical cache

immutable OpticalCache{T,NDri,NDec,Tra}
    drives::NTuple{NDri,TrigCache{T}}
    decays::NTuple{NDec,TrigCache{T}}
    trackers::Tra
end

@compat @inline function (::Type{OpticalCache{T,NDri,NDec}}){T,NDri,NDec,Tra}(drives, decays, trackers::Tra)
    OpticalCache{T,NDri,NDec,Tra}(drives, decays, trackers)
end

@generated function OpticalCache{M<:MotionSystem}(sys::M, _dx, _dt, nele)
    T = System.get_value_type(M)
    NDri = length(System.get_drive_types(M))
    NDec = length(System.get_transition_types(M))

    init_expr = quote
        dx = $T(_dx)
        dt = $T(_dt)
        nele_2 = nele ÷ 2
        xs = ((1:nele) - nele_2) * dx
    end

    calc_drives = quote
        drives = ($([:(TrigCache(sys.drives[$i], xs)) for i in 1:NDri]...),)
    end
    calc_decays = quote
        decays = ($([:(TrigCache(sys.intern.transitions[$i], xs))
                     for i in 1:NDec]...),)
    end
    calc_trackers = quote
        trackers = ($([:(Optical.PhaseTracker(sys.drives[$i]))
                       for i in 1:NDri]...),)
    end

    quote
        $init_expr
        $calc_drives
        $calc_decays
        $calc_trackers
        OpticalCache{$T,$NDri,$NDec}(drives, decays, trackers)
    end
end

# Coupling cache

"""
Discribes the strength of a drive on a certain transition.
Determined by the drive, the transition and the quantization axis
"""
immutable Coupling{T}
    P_off::Complex{T} # i * e^(iϕ) * sin(Ω dt)
    P_diag::T # cos(Ω dt)
    # Where ϕ is the initial phase of the coupling, Ω is the Rabi frequency
end

immutable CouplingCache{T,N,Idxs}
    # Idxs:
    #     Tuple of (drive_id, transition_id)
    couplings::NTuple{N,Coupling{T}}
end

immutable DriveCoupling{T}
    amp::T
    overlap_σ⁺::Complex{T}
    overlap_σ⁻::Complex{T}
    overlap_π::Complex{T}

    couple_σ⁺::Bool
    couple_σ⁻::Bool
    couple_π::Bool

    function DriveCoupling(ax::Vec3D{T}, _amp::Vec3D)
        abs_amp = abs(_amp)
        amp = _amp / abs_amp
        overlap_σ⁺ = (ax, Trans_σ⁺) * amp
        overlap_σ⁻ = (ax, Trans_σ⁻) * amp
        overlap_π = (ax, Trans_π) * amp

        new(abs_amp, overlap_σ⁺, overlap_σ⁻, overlap_π,
            abs(overlap_σ⁺) >= 1e-3, abs(overlap_σ⁻) >= 1e-3,
            abs(overlap_π) >= 1e-3)
    end
end

Base.getindex(cpl::DriveCoupling, trans::TransitionType) = if trans == Trans_σ⁺
    (cpl.couple_σ⁺, cpl.overlap_σ⁺)
elseif trans == Trans_σ⁻
    (cpl.couple_σ⁻, cpl.overlap_σ⁻)
else
    (cpl.couple_π, cpl.overlap_π)
end

@generated function CouplingCache{M<:MotionSystem}(sys::M, _dx, _dt, nele)
    T = System.get_value_type(M)
    idxs = NTuple{2,Int}[]
    _couplings = Expr[]

    init_expr = quote
        dt = $T(_dt)
        intern = sys.intern
        transitions = intern.transitions
    end

    Tdrives = System.get_drive_types(M)
    Ttrans = System.get_transition_types(M)
    trans_pairs = System.get_transition_pairs(M)
    Ax = System.get_quant_axis(M)

    drive_gids = System.get_drive_gids(M)
    state_gids = System.get_state_gids(M)
    # snames = System.get_state_names(M)

    for i in 1:length(Tdrives)
        drive_from, drive_to = drive_gids[i]
        Tdrive = Tdrives[i]
        drive_amp = Optical.get_drive_amplitude(Tdrive)
        cpl = DriveCoupling{T}(Ax, drive_amp)
        for j in 1:length(Ttrans)
            from_id, to_id = trans_pairs[j]
            from_grp = state_gids[from_id]
            to_grp = state_gids[to_id]
            ((from_grp, to_grp) == (drive_from, drive_to) ||
             (from_grp, to_grp) == (drive_to, drive_from)) || continue

            Ttran = Ttrans[j]
            has_couple, overlap = cpl[Atomic.get_transition_type(Ttran)]
            has_couple || continue
            amp_eff = overlap * cpl.amp
            push!(idxs, (i, j))

            # println((i, j, snames[from_id], snames[to_id]))

            ex = quote
                Ω = $amp_eff * transitions[$j].α
                expϕ0 = sign(Ω)
                dθ = abs(Ω) * dt
                sindθ = sin(dθ)
                cosdθ = cos(dθ)
                Coupling{$T}(((1im) * expϕ0) * sindθ, cosdθ)
            end
            push!(_couplings, Expr(:let, ex))
        end
    end
    # println()

    N = length(idxs)
    Idxs = (idxs...)
    quote
        $init_expr
        CouplingCache{$T,$N,$Idxs}(($(_couplings...),))
    end
end

immutable SystemPropagator{Sys,T,Mc,Oc,Cc,P,PI}
    sys::Sys

    dt::T
    dx::T
    nstep::Int
    nele::Int

    motion::Mc
    optical::Oc
    coupling::Cc

    tmp::Matrix{Complex{T}} # nele x nstates
    sotmp::SoCMatrix{T} # nele x nstates

    p_fft!::P
    p_bfft!::PI
end

function SystemPropagator{Sys<:MotionSystem}(sys::Sys, _dt, _dx, nstep, nele)
    T = System.get_value_type(Sys)
    dt = T(_dt)
    dx = T(_dx)

    motion = HMotionCache(sys, dx, dt, nele)
    Mc = typeof(motion)
    optical = OpticalCache(sys, dx, dt, nele)
    Oc = typeof(optical)
    coupling = CouplingCache(sys, dx, dt, nele)
    Cc = typeof(coupling)

    nstates = System.num_states(Sys)
    tmp = Matrix{Complex{T}}(nele, nstates)
    sotmp = convert(StructOfArrays, tmp)

    p_fft! = plan_fft!(tmp, 1, flags=FFTW.MEASURE)
    p_bfft! = plan_bfft!(tmp, 1, flags=FFTW.MEASURE)
    P = typeof(p_fft!)
    PI = typeof(p_bfft!)
    SystemPropagator{Sys,T,Mc,Oc,Cc,P,PI}(sys, dt, dx, nstep, nele, motion,
                                          optical, coupling, tmp, sotmp,
                                          p_fft!, p_bfft!)
end

include("propagate.jl")

end
