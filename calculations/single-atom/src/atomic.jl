#!/usr/bin/julia -f

###
# Internal states of the atom

@enum Polarization Pol_σplus Pol_σminus Pol_π

immutable Amplitude3D{T}
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

@inline call{T}(::Type{Amplitude3D}, pol::Polarization, amp::T) =
    Amplitude3D{T}(pol, amp)
@inline call{T}(::Type{Amplitude3D}, pol::Polarization, amp::Complex{T}) =
    Amplitude3D{T}(pol, amp)
