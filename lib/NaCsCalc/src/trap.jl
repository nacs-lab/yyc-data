#!/usr/bin/julia

module Trap

import GSL

function η(m, freq, k)
    ħ = 1.0545718e-34
    z_0 = √(ħ / 2 / m / (2π * freq))
    z_0 * k
end

function sideband{T1,T2,T3}(n1::T1, n2::T2,
                            η::T3)::float(promote_type(T1, T2, T3))
    if n1 < 0 || n2 < 0
        return η
    elseif η == 0
        if n1 == n2
            return 1
        else
            return 0
        end
    end
    # Ref http://journals.aps.org/pra/pdf/10.1103/PhysRevA.20.1521
    # Δn ≡ |n1 - n2|
    # n₋ ≡ min(n1, n2)
    # n₊ ≡ max(n1, n2)
    #   ⟨n1|exp(ikx)|n2⟩
    # = ⟨n1|exp(iη(a + a†))|n2⟩
    # = exp(-η^2 / 2) η^Δn √(γ(n₋ + 1) / γ(n₊ + 1)) L^Δn_n₋(η^2)
    # = exp(-η^2 / 2 + Δn log(η) + lγ(n₋ + 1) / 2 - lγ(n₊ + 1) / 2) L^Δn_n₋(η^2)
    n₋ = min(n1, n2)
    n₊ = max(n1, n2)
    Δn = abs(n1 - n2)
    η² = η * η
    lpre = (-η² + lgamma(n₋ + 1) - lgamma(n₊ + 1)) / 2 + log(η) * Δn
    lag = GSL.sf_laguerre_n(n₋, Δn, η²)
    lag * exp(lpre)
end

end
