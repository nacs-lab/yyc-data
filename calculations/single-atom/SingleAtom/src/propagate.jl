#!/usr/bin/julia -f

module Propagate

using ..Utils
using ..System
using ..Atomic
import ..Optical

immutable HMotionCache{T,N}
    E_k::Vector{T}
    P_k::SoCVector{T}

    E_x::NTuple{N,Vector{T}}
    P_x2::NTuple{N,SoCVector{T}}
end

@generated function HMotionCache{M<:MotionSystem}(sys::M, _dx, _dt, nele)
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
            P_k[i] = exp(im * e_k * dt)
        end
    end
    e_temps = [gensym() for i in 1:NPots]
    calc_x = quote
        @inbounds for i in 1:nele
            x = (i - nele_2) * dx
            $([:($(e_temps[i]) = get_potential(sys.potentials[$i], m, x))
               for i in 1:NPots]...)
            $([:(E_x[$i][i] = $(e_temps[i])) for i in 1:NPots]...)
            $([:(P_x2[$i][i] = exp(im * dt / 2 * $(e_temps[i])))
               for i in 1:NPots]...)
        end
    end
    quote
        $init_expr
        $calc_k
        $calc_x
        HMotionCache{$T,$NPots}(E_k, P_k, E_x, P_x2)
    end
end

# Optical cache

immutable OpticalCache{T,NDri,NDec,Tra}
    drives::NTuple{NDri,TrigCache{T}}
    decays::NTuple{NDec,TrigCache{T}}
    trackers::Tra
end

@inline function call{T,NDri,NDec,Tra}(::Type{OpticalCache{T,NDri,NDec}},
                                       drives, decays, trackers::Tra)
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
    Ω::T # Rabi rate
    sindθ::T # sin(Ω dt)
    cosdθ::T # cos(Ω dt)
end

immutable CouplingCache{T,N,Idxs}
    # Idxs:
    #     Tuple of (drive_id, transition_id)
    couplings::NTuple{N,Coupling{T}}
end

immutable DriveCoupling{T}
    amp::T
    overlap_σ⁺::T
    overlap_σ⁻::T
    overlap_π::T

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
            overlap_σ⁺ >= 1e-3, overlap_σ⁻ >= 1e-3, overlap_π >= 1e-3)
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
    Ax = System.get_quant_axis(M)

    for i in 1:length(Tdrives)
        Tdrive = Tdrives[i]
        drive_amp = Optical.get_drive_amplitude(Tdrive)
        cpl = DriveCoupling{T}(Ax, drive_amp)
        for j in 1:length(Ttrans)
            Ttran = Ttrans[j]
            has_couple, overlap = cpl[Atomic.get_transition_type(Ttran)]
            has_couple || continue
            amp_eff = overlap * cpl.amp
            push!(idxs, (i, j))
            ex = quote
                Ω = $amp_eff * transitions[$j].α
                dθ = Ω * dt
                sindθ = sin(dθ)
                cosdθ = cos(dθ)
                Coupling{$T}(Ω, sindθ, cosdθ)
            end
            push!(_couplings, Expr(:let, ex))
        end
    end

    N = length(idxs)
    Idxs = (idxs...)
    quote
        $init_expr
        CouplingCache{$T,$N,$Idxs}(($(_couplings...),))
    end
end

end
