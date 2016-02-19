#!/usr/bin/julia -f

module Sidebands

using GSL

function strength(n1::Integer, n2::Integer, η)
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
    lag = sf_laguerre_n(n₋, Δn, η²)
    lag * exp(lpre)
end

function thermal_strength(nbar, Δn, η)
    s = 0.0
    nmin = max(0, -Δn)
    nmax = nmin + round(Int, nbar * 10)
    weight = 0.0
    for n in 0:nmax
        w = exp(-n / nbar)
        weight += w
        n < nmin && continue
        s += strength(n, n + Δn, η)^2 * w
    end
    s / weight
end

end
