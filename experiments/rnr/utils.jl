#

using SpecialFunctions

function G(x::Float64, t::Float64)
    x1 = x * (2 + t^2)
    x2 = x * t * √(4 + t^2)
    @fastmath if x1 <= 400 && x2 <= 400
        return 2 * exp(-x1) * besseli(0, x2)
    else
        z8 = 1 / 8 / x2
        expansion = @evalpoly(z8, 1, 1, 9 / 2, 9 * 25 / 6, 9 * 25 * 49 / 24)
        return 2 * exp(x2 - x1) / √(2π * x2) * expansion
    end
end
G(x, t) = G(Float64(x), Float64(t))

_tuple_tail(x, y...) = y

# arguments: k_BT, time, energy, trapping frequency
#     units: T^-1,    T,   T^-1,               T^-1
@inline function pdf_rnr_1d(kT, t0, ω, E)
    β = ω / kT
    ɛ = E / ω
    t = ω * t0
    tf = @fastmath tanh(β / 2)
    return tf * G(ɛ * tf, t) / ω
end

@inline function pdf_rnr{N}(kT, t, ωs::NTuple{N}, Es::NTuple{N})
    return pdf_rnr_1d(kT, t, ωs[1], Es[1]) * pdf_rnr(kT, t, _tuple_tail(ωs...),
                                                      _tuple_tail(Es...))
end

@inline function pdf_rnr(kT, t, ωs::NTuple{1}, Es::NTuple{1})
    return pdf_rnr_1d(kT, t, ωs[1], Es[1])
end

function pdf_rnr_polar(kT, t, ωs, Er, θ, ϕ)
    sinθ = @fastmath sin(θ)
    cosθ = @fastmath cos(θ)
    sinϕ = @fastmath sin(ϕ)
    cosϕ = @fastmath cos(ϕ)

    Ex = Er * cosθ
    Eρ = Er * sinθ
    Ey = Eρ * cosϕ
    Ez = Eρ * sinϕ
    return Er^2 * pdf_rnr(kT, t, ωs, (Ex, Ey, Ez)) * sinθ
end
