#!/usr/bin/julia

import GSL

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

coupling_D(D2::Bool, Ix2, F1x2, F2x2, mF1x2, mF2x2) =
    coupling_LS(0, 2, 1, D2 ? 3 : 1, Ix2, F1x2, F2x2, mF1x2, mF2x2)

function na_scattering_d2(F1::Integer, mF1::Integer, pol::Integer, F2::Integer, mF2::Integer)
    if !(F1 in (1, 2)) || !(F2 in (1, 2))
        throw(ArgumentError("Invalid F for Na"))
    elseif !(pol in (-1, 0, 1))
        throw(ArgumentError("Invalid polarization"))
    end
    mF′ = mF1 + pol

    F1x2 = F1 * 2
    mF1x2 = mF1 * 2
    F2x2 = F2 * 2
    mF2x2 = mF2 * 2
    mF′x2 = mF′ * 2

    cpl = 0.0
    for F′ in (F1 - 1):(F1 + 1)
        if abs(mF′) > F′
            continue
        end
        if abs(F′ - F2) > 1
            continue
        end
        F′x2 = F′ * 2
        cpl += (coupling_D(true, 3, F1x2, F′x2, mF1x2, mF′x2) *
                coupling_D(true, 3, F2x2, F′x2, mF2x2, mF′x2))
    end
    return cpl
end
