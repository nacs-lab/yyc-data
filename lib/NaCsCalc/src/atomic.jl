#!/usr/bin/julia

module Atomic

import GSL

"""
    coupling_LS(L1x2, L2x2, J1x2, J2x2, Ix2, F1x2, F2x2, mF1x2, mF2x2)

Compute the relative matrix element between two states (`F1`, `mF1`) and (`F2`, `mF2`)
relative to the reduced matrix element on two `L` states for a single electron atom
in the LS coupling regime (i.e. Hydrogen or Alkali-like).
The `I` for the atom and `L`, `J`, `F` and `mF` for the two states are given in their value
multiplying by two so they are all integers.
"""
function coupling_LS(L1x2::Integer, L2x2::Integer, J1x2::Integer, J2x2::Integer, Ix2::Integer,
                     F1x2::Integer, F2x2::Integer, mF1x2::Integer, mF2x2::Integer)::Float64
    if abs(L1x2 - L2x2) != 2
        return 0
    elseif abs(J1x2 - L1x2) != 1 || abs(J2x2 - L2x2) != 1
        throw(ArgumentError("Illegal J"))
    elseif !(abs(J1x2 - Ix2) <= F1x2 <= J1x2 + Ix2 && abs(J2x2 - Ix2) <= F2x2 <= J2x2 + Ix2)
        throw(ArgumentError("Illegal F"))
    elseif (F1x2 + J1x2 + Ix2) % 2 != 0 || (F2x2 + J2x2 + Ix2) % 2 != 0
        throw(ArgumentError("Illegal F"))
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
    rates = Matrix{Float64}(nstates, nstates)
    @inbounds for i in 1:nstates
        r, F1x2, mF1x2 = idx_to_state(i)
        for j in 1:nstates
            _, F2x2, mF2x2 = idx_to_state(j)
            rates[j, i] = scatter_D(D2, Ix2, F1x2, F2x2, mF1x2, mF2x2, pol) * r
        end
    end
    return rates
end

function all_scatter_D(D2::Bool, Ix2, pol, fscale)
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
    rates = Matrix{Float64}(nstates, nstates)
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

end
