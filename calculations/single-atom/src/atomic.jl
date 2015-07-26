#!/usr/bin/julia -f

###
# Internal states of the atom

@enum Polarization Pol_σplus Pol_σminus Pol_π

abstract Amplitude

Base.abs(amp::Amplitude) = sqrt(abs2(amp))

immutable Amplitude3D{T<:Real} <: Amplitude
    x::Complex{T}
    y::Complex{T}
    z::Complex{T}
    @inline Amplitude3D(x, y, z) = new(x, y, z)
    @inline function Amplitude3D(pol::Polarization,
                                 amp::Union{T, Complex{T}}=T(1))
        if pol == Pol_π
            return new(0, 0, amp)
        end
        amp2 = sqrt(T(0.5)) * amp
        if pol == Pol_σplus
            new(amp2, im * amp2, 0)
        else
            new(amp2, -im * amp2, 0)
        end
    end
end

@inline call{T<:Real}(::Type{Amplitude3D}, pol::Polarization, amp::T) =
    Amplitude3D{T}(pol, amp)

@inline call{T<:Real}(::Type{Amplitude3D}, pol::Polarization,
                      amp::Complex{T}) = Amplitude3D{T}(pol, amp)

Base.abs2(amp::Amplitude3D) = abs2(amp.x) + abs2(amp.y) + abs2(amp.z)

immutable LinearPol
end
immutable CircularPol{Dir}
end

const linearPol = LinearPol()
const rightCirPol = CircularPol{1}()
const leftCirPol = CircularPol{-1}()

immutable Amplitude2D{T<:Real} <: Amplitude
    x::Complex{T}
    y::Complex{T}
    @inline Amplitude2D(x, y) = new(x, y)
    @inline Amplitude2D(::LinearPol, amp, angl=T(0)) =
        new(amp * cos(angl), amp * sin(angl))
    @inline function Amplitude2D{Dir}(::CircularPol{Dir}, amp)
        amp2 = sqrt(T(0.5)) * amp
        new(amp2, Dir * im * amp2)
    end
end

@inline call{T}(::Type{Amplitude2D}, p::LinearPol, amp::T, angl=T(0)) =
    Amplitude2D{T}(p, amp, angl)

@inline call{T}(::Type{Amplitude2D}, p::LinearPol, amp::Complex{T},
                angl=T(0)) = Amplitude2D{T}(p, amp, angl)

@inline call{T}(::Type{Amplitude2D}, p::CircularPol, amp::T) =
    Amplitude2D{T}(p, amp)

@inline call{T}(::Type{Amplitude2D}, p::CircularPol, amp::Complex{T}) =
    Amplitude2D{T}(p, amp)

Base.abs2(amp::Amplitude2D) = abs2(amp.x) + abs2(amp.y)
