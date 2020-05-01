#!/usr/bin/julia

module Atomic

import GSL
import ..Utils

"""
    coupling_LS(L1x2, L2x2, J1x2, J2x2, Ix2, F1x2, F2x2, mF1x2, mF2x2)

Compute the relative matrix element between two states (`F1`, `mF1`) and (`F2`, `mF2`)
relative to the reduced matrix element on two `L` states for a single electron atom
in the LS coupling regime (i.e. Hydrogen or Alkali-like).
The `I` for the atom and `L`, `J`, `F` and `mF` for the two states are given in their value
multiplying by two so they are all integers.

Note that the reduced matrix element between L=0 and L=1 is `√(3)` that of the
actual matrix element for each transition. This is a result of following the steck
convention.
"""
function coupling_LS(L1x2::Integer, L2x2::Integer, J1x2::Integer, J2x2::Integer, Ix2::Integer,
                     F1x2::Integer, F2x2::Integer, mF1x2::Integer, mF2x2::Integer)::Float64
    if abs(L1x2 - L2x2) != 2
        return 0
    elseif abs(J1x2 - L1x2) != 1 || abs(J2x2 - L2x2) != 1
        throw(ArgumentError("Illegal J"))
    elseif (F1x2 + J1x2 + Ix2) % 2 != 0 || (F2x2 + J2x2 + Ix2) % 2 != 0
        throw(ArgumentError("Illegal F"))
    elseif !(abs(J1x2 - Ix2) <= F1x2 <= J1x2 + Ix2 && abs(J2x2 - Ix2) <= F2x2 <= J2x2 + Ix2)
        return 0
    elseif !((F1x2 - F2x2) in (-2, 0, 2) && (mF1x2 - mF2x2) in (-2, 0, 2))
        return 0
    end
    qx2 = mF1x2 - mF2x2
    m1powx2 = F2x2 + mF1x2 - 2
    factor = sqrt(F1x2 + 1) * GSL.sf_coupling_3j(F2x2, 2, F1x2, mF2x2, qx2, -mF1x2)
    m1powx2 += F2x2 + J1x2 + 2 + Ix2
    factor *= sqrt((F2x2 + 1) * (J1x2 + 1)) * GSL.sf_coupling_6j(J1x2, J2x2, 2, F2x2, F1x2, Ix2)
    m1powx2 += J2x2 + L1x2 + 2 + 1
    factor *= sqrt((J2x2 + 1) * (L1x2 + 1)) * GSL.sf_coupling_6j(L1x2, L2x2, 2, J2x2, J1x2, 1)
    if m1powx2 % 4 == 2
        factor = -factor
    end
    return factor
end

"""
    coupling_D(D2::Bool, Ix2, F1x2, F2x2, mF1x2, mF2x2)

Compute the relative matrix element between two states (`F1`, `mF1`) and (`F2`, `mF2`)
on the D line of a Alkali-like atom, i.e. `s` to `p` transition.
"""
coupling_D(D2::Bool, Ix2, F1x2, F2x2, mF1x2, mF2x2) =
    coupling_LS(0, 2, 1, D2 ? 3 : 1, Ix2, F1x2, F2x2, mF1x2, mF2x2)

@inline dummy_fscale(F, mF) = 1

"""
    coupling_D_raman(D2::Bool, Ix2, F1x2, F2x2, mF1x2, mF2x2, pol)

Compute the Raman coupling between two Zeeman ground state with off-resonance Raman transition
on a D line. `pol ∈ (-1, 0, 1)` specifies the polarization of the beam addressing the first state.
"""
function coupling_D_raman(D2::Bool, Ix2::Integer, F1x2::Integer, F2x2::Integer,
                          mF1x2::Integer, mF2x2::Integer, pol::Integer, fscale=dummy_fscale)
    pol in (-1, 0, 1) || throw(ArgumentError("Invalid polarization"))
    mF′x2 = mF1x2 + pol * 2

    cpl = 0.0
    for F′x2 in (F1x2 - 2):2:(F1x2 + 2)
        if abs(mF′x2) > F′x2
            continue
        end
        if abs(F′x2 - F2x2) > 2
            continue
        end
        cpl += (coupling_D(D2, Ix2, F1x2, F′x2, mF1x2, mF′x2) *
                coupling_D(D2, Ix2, F2x2, F′x2, mF2x2, mF′x2)) * fscale(F′x2, mF′x2)
    end
    return cpl
end

function scatter_bothD(Ix2::Integer, F1x2::Integer, F2x2::Integer,
                       mF1x2::Integer, mF2x2::Integer, pol::NTuple{3,Real},
                       fscale=(D2, F, mF)->1)
    rate = 0.0
    for i in 1:3
        p = pol[i]
        if p == 0
            continue
        end
        rate += abs2(coupling_D_raman(false, Ix2, F1x2, F2x2, mF1x2, mF2x2, i - 2,
                                      (F, mF)->fscale(false, F, mF)) +
                     coupling_D_raman(true, Ix2, F1x2, F2x2, mF1x2, mF2x2, i - 2,
                                      (F, mF)->fscale(true, F, mF))) * p
    end
    return rate
end

function scatter_D(D2::Bool, Ix2::Integer, F1x2::Integer, F2x2::Integer,
                   mF1x2::Integer, mF2x2::Integer, pol::NTuple{3,Real}, fscale=dummy_fscale)
    rate = 0.0
    for i in 1:3
        p = pol[i]
        if p == 0
            continue
        end
        rate += abs2(coupling_D_raman(D2, Ix2, F1x2, F2x2, mF1x2, mF2x2, i - 2, fscale)) * p
    end
    return rate
end

function all_scatter_D(D2::Bool, Ix2, pol, rhi, rlo)
    # Order: F high->low; mF low->high
    idx_to_state = function (idx)
        local Fx2, mFx2, r
        if idx > Ix2 + 2
            # low F
            idx -= Ix2 + 2
            Fx2 = Ix2 - 1
            mFx2 = idx * 2 - Ix2 - 1
            r = rlo
        else
            # high F
            Fx2 = Ix2 + 1
            mFx2 = idx * 2 - Ix2 - 3
            r = rhi
        end
        return r, Fx2, mFx2
    end
    nstates = (Ix2 + 1) * 2
    rates = Matrix{Float64}(undef, nstates, nstates)
    @inbounds for i in 1:nstates
        r, F1x2, mF1x2 = idx_to_state(i)
        for j in 1:nstates
            _, F2x2, mF2x2 = idx_to_state(j)
            rates[j, i] = scatter_D(D2, Ix2, F1x2, F2x2, mF1x2, mF2x2, pol) * r
        end
    end
    return rates
end

function all_scatter_D(fscale, D2::Bool, Ix2, pol)
    # Order: F high->low; mF low->high
    idx_to_state = function (idx)
        local Fx2, mFx2
        if idx > Ix2 + 2
            # low F
            idx -= Ix2 + 2
            Fx2 = Ix2 - 1
            mFx2 = idx * 2 - Ix2 - 1
        else
            # high F
            Fx2 = Ix2 + 1
            mFx2 = idx * 2 - Ix2 - 3
        end
        return Fx2, mFx2
    end
    nstates = (Ix2 + 1) * 2
    rates = Matrix{Float64}(undef, nstates, nstates)
    @inbounds for i in 1:nstates
        F1x2, mF1x2 = idx_to_state(i)
        for j in 1:nstates
            let (F2x2, mF2x2) = idx_to_state(j)
                rates[j, i] = scatter_D(D2, Ix2, F1x2, F2x2, mF1x2, mF2x2, pol,
                                        (F, mF)->fscale(F1x2, mF1x2, F, mF))
            end
        end
    end
    return rates
end

function all_scatter_bothD(fscale, Ix2, pol)
    # Order: F high->low; mF low->high
    idx_to_state = function (idx)
        local Fx2, mFx2
        if idx > Ix2 + 2
            # low F
            idx -= Ix2 + 2
            Fx2 = Ix2 - 1
            mFx2 = idx * 2 - Ix2 - 1
        else
            # high F
            Fx2 = Ix2 + 1
            mFx2 = idx * 2 - Ix2 - 3
        end
        return Fx2, mFx2
    end
    nstates = (Ix2 + 1) * 2
    rates = Matrix{Float64}(undef, nstates, nstates)
    @inbounds for i in 1:nstates
        F1x2, mF1x2 = idx_to_state(i)
        for j in 1:nstates
            let (F2x2, mF2x2) = idx_to_state(j)
                rates[j, i] = scatter_bothD(Ix2, F1x2, F2x2, mF1x2, mF2x2, pol,
                                            (D2, F, mF)->fscale(D2, F1x2, mF1x2, F, mF))
            end
        end
    end
    return rates
end

struct Alkali
    Ix2::Int
    D1::Float64
    D2::Float64
    Γ1::Float64
    Γ2::Float64
    HF_S::NTuple{2,Float64}
    HF_P1_2::NTuple{2,Float64}
    HF_P3_2::NTuple{4,Float64}
end

function _check_j1_2(atom::Alkali, name, Fx2)
    if Fx2 != atom.Ix2 + 1 && Fx2 != atom.Ix2 - 1
        throw(ArgumentError(
            "Invalid F for $(name)1/2 state with I=$(atom.Ix2 / 2): $(Fx2 / 2)"))
    end
end

function _check_j3_2(atom::Alkali, Fx2)
    if Fx2 != atom.Ix2 + 1 && Fx2 != atom.Ix2 - 1 && Fx2 != atom.Ix2 + 3 && Fx2 != atom.Ix2 - 3
        throw(ArgumentError(
            "Invalid F for P3/2 state with I=$(atom.Ix2 / 2): $(Fx2 / 2)"))
    end
end

function _check_state(atom::Alkali, Lx2, Jx2, Fx2, mFx2=Fx2)
    if Lx2 < 0 || Jx2 < 0 || Fx2 < 0
        throw(ArgumentError("Negative angular momentum"))
    elseif Lx2 == 0
        Jx2 == 1 || throw(ArgumentError("Invalid J for S state: $(Jx2 / 2)"))
        _check_j1_2(atom, :S, Fx2)
    elseif Lx2 == 2
        if Jx2 == 1
            _check_j1_2(atom, :P, Fx2)
        elseif Jx2 == 3
            _check_j3_2(atom, Fx2)
        else
            throw(ArgumentError("Invalid J for P state: $(Jx2 / 2)"))
        end
    else
        throw(ArgumentError("Invalid L $(Lx2 / 2)"))
    end
    if mFx2 > Fx2 || mFx2 < -Fx2 || (Fx2 - mFx2) % 2 != 0
        throw(ArgumentError("Invalid m_F for F=$(Fx2 / 2): $(mFx2 / 2)"))
    end
end

function _get_freq_s(atom, Jx2, Fx2)
    Fx2 == atom.Ix2 + 1 && return atom.HF_S[1]
    return atom.HF_S[2]
end

function _get_freq_d1(atom, Fx2)
    Fx2 == atom.Ix2 + 1 && return atom.D1 + atom.HF_P1_2[1]
    return atom.D1 + atom.HF_P1_2[2]
end

function _get_freq_d2(atom, Fx2)
    Fx2 == atom.Ix2 + 3 && return atom.D2 + atom.HF_P3_2[1]
    Fx2 == atom.Ix2 + 1 && return atom.D2 + atom.HF_P3_2[2]
    Fx2 == atom.Ix2 - 1 && return atom.D2 + atom.HF_P3_2[3]
    return atom.D2 + atom.HF_P3_2[4]
end

function _get_freq_p(atom, Jx2, Fx2)
    Jx2 == 1 && return _get_freq_d1(atom, Fx2)
    Jx2 == 3 && return _get_freq_d2(atom, Fx2)
end

function _get_freq(atom::Alkali, Lx2, Jx2, Fx2)
    Lx2 == 0 && return _get_freq_s(atom, Jx2, Fx2)
    return _get_freq_p(atom, Jx2, Fx2)
end

@inline function get_freq(atom::Alkali, Lx2, Jx2, Fx2)
    @boundscheck _check_state(atom, Lx2, Jx2, Fx2)
    _get_freq(atom, Lx2, Jx2, Fx2)
end

function _get_amp_j(atom, F1x2, mF1x2, F2x2, mF2x2, Jex2, freq0, p, sqrtΓ)::Float64
    get_scale = function (Fex2, mFex2)
        @inbounds freqe = get_freq(atom, 2, Jex2, Fex2)
        return sqrtΓ / (2π) / (freqe - freq0)
    end
    return coupling_D_raman(Jex2 == 3, atom.Ix2, F1x2, F2x2, mF1x2, mF2x2, p, get_scale)
end

function _get_amp_p(atom, F1x2, mF1x2, F2x2, mF2x2, freq0, p)::Float64
    return _get_amp_j(atom, F1x2, mF1x2, F2x2, mF2x2, 1, freq0, p, sqrt(atom.Γ1)) +
        _get_amp_j(atom, F1x2, mF1x2, F2x2, mF2x2, 3, freq0, p, sqrt(atom.Γ2))
end

# This function include the effect of coherence between all excited states in P1/2 and P3/2.
# This only work in low power limit where the excited state population can be ignored
# (the transition is not saturated or the scattering rate is much lower than the natural decay
# rate of the transition)
function get_scatter(atom::Alkali, Ω, freq, F1x2, mF1x2, F2x2, mF2x2, pol)
    _check_state(atom, 0, 1, F1x2, mF1x2)
    _check_state(atom, 0, 1, F2x2, mF2x2)
    @inbounds freqg = get_freq(atom, 0, 1, F1x2)
    freq0 = freqg + freq
    rate = 0.0
    for i in 1:3
        p = pol[i]
        if p == 0
            continue
        end
        rate += abs2(_get_amp_p(atom, F1x2, mF1x2, F2x2, mF2x2, freq0, i - 2)) * p
    end
    # The 9 below is the factor of 3 (√3 in each leg) in amplitude coming from the
    # reduced matrix element of L above.
    return rate * Ω^2 / 8 * 9
end

end
