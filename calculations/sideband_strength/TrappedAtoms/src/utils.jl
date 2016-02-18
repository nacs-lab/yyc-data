#!/usr/bin/julia -f

function calculate_η(m, freq, k)
    ħ = 1.0545718e-34
    z_0 = √(ħ / 2 / m / (2π * freq))
    z_0 * k
end
