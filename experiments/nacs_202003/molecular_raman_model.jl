#!/usr/bin/julia

function evolve(t, δ, Ω, Γ1, Γ2)
    M = [Γ1 + im * δ Ω
         -Ω Γ2 - im * δ]
    M .= M .* (-t/2)
    return abs2((exp(M) * [1, 0])[1]) # expm
end

function model_2d(t, f, p)
    p0, p1, f0, Ω, Γ1, Γ2 = p
    return p0 + p1 * evolve(t, 2π * (f - f0), Ω, Γ1, Γ2)
end
