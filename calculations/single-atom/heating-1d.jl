#!/usr/bin/julia -f

# 1D harmonic trap with magic wavelength using Monte Carlo

immutable MagicHarmonic1D{T}
    ω::T
    center::T

    k_emit::T
    Γ::T

    k_drive::T
    Ω::T
    δ::T
end

immutable Propagator1D{T}
    H::MagicHarmonic1D{T}
    dt::T
    dx::T
    nstep::Int
    nele::Int

    P_k::Vector{Complex{T}} # nele
    P_x2::Vector{Complex{T}} # nele
    P_σ11::Matrix{Complex{T}} # nele x nstep
    P_σ12::Matrix{Complex{T}} # nele x nstep
    # P_σ21::Matrix{Complex{T}} # nele x nstep = H_σ12*
    P_σ22::Matrix{Complex{T}} # nele x nstep

    function Propagator1D(H::MagicHarmonic1D{T}, dt, dx, nstep, nele)
        # Hamiltonian is
        #   H = k^2 / 2 + ω^2 x^2 / 2
        #       + Ω (cos(δ t + k_drive x) σ_x + sin(δ t + k_drive x) σ_y)
        #     = H_p + H_x + H_σ
        # H_k = k^2 / 2
        # H_x = ω^2 x^2 / 2
        # H_σ = Ω (cos(δ t + k_drive x) σ_x + sin(δ t + k_drive x) σ_y)

        # Propagator is
        # P_k = exp(im H_k Δt) = exp(im k^2 Δt / 2)
        # P_x = exp(im H_x Δt) = exp(im ω^2 x^2 Δt / 2)
        # P_k2 = exp(im H_x Δt / 2) = exp(im ω^2 x^2 Δt / 4)

        P_k = Vector{Complex{T}}(nele)
        P_x2 = Vector{Complex{T}}(nele)

        k0 = 2π / (nele * dx)
        @inbounds for i in 1:nele
            x = i * dx # coordinate
            P_x2[i] = exp(im * H.ω^2 * x^2 * dt / 4)

            k1 = i - 1
            k2 = k1 - nele
            k = ifelse(k1 + k2 <= 0, k1, k2) * k0
            P_k[i] = exp(im * k^2 * dt / 2)
        end

        # P_σ = exp(im H_σ Δt)
        #     = exp(im Ω (cos(δ t + k_drive x) σ_x +
        #                 sin(δ t + k_drive x) σ_y) Δt)
        #     = cos(Ω Δt) + im * (cos(δ t + k_drive x) σ_x +
        #                         sin(δ t + k_drive x) σ_y) * sin(Ω Δt)
        #     = cos(Ω Δt) + im cos(δ t + k_drive x) sin(Ω Δt) σ_x +
        #       im sin(δ t + k_drive x) sin(Ω Δt) σ_y
        #     = [cos(Ω Δt), im exp(im(δ t + k_drive x)) sin(Ω Δt)
        #        im exp(-im(δ t + k_drive x)) sin(Ω Δt), cos(Ω Δt)]

        P_σ11 = Matrix{Complex{T}}(nele, nstep)
        P_σ12 = Matrix{Complex{T}}(nele, nstep)
        P_σ22 = Matrix{Complex{T}}(nele, nstep)

        cos_dt = cos(Ω * dt)
        sin_dt = sin(Ω * dt)
        @inbounds for j in 1:nstep
            t = i * dt
            θ_t = H.δ * t
            for i in 1:nele
                x = i * dx # coordinate
                θ_x = H.k_drive * x
                θ = θ_t + θ_x

                # P_σ12[i, j] = im exp(im θ) sin_dt
                #             = (im cos(θ) - sin(θ)) sin_dt

                P_σ11[i, j] = P_σ22[i, j] = cos_dt
                P_σ12[i, j] = (im * cos(θ) - sin(θ)) * sin_dt
            end
        end

        new(H, dt, dx, nstep, nele, P_k, P_x2, P_σ11, P_σ12, P_σ22)
    end
end

function propagate{T}(P::Propagator1D{T},
                      ψ0::Matrix{Complex{T}}, # 2 x nele
                      ψs::Array{Complex{T}, 3} # 2 x nele x (nstep + 1)
                      )
    @inbounds ψs[:, :, 1] = ψ0

    tmp = similar(ψ0)

    # FFT plan
    p_fft! = plan_fft!(tmp, 2, FFTW.MEASURE)
    p_ifft! = plan_ifft!(tmp, 2, FFTW.MEASURE)

    @inbounds for i in 2:nstep
        for j in 1:nele
            p_x2 = P.P_x2[j]
            tmp[1, j] = ψs[1, j, i - 1] * p_x2
            tmp[2, j] = ψs[2, j, i - 1] * p_x2
        end
        p_fft!(tmp)
        for j in 1:nele
            p_k = P.P_k[j]
            tmp[1, j] *= p_k
            tmp[2, j] *= p_k
        end
        p_ifft!(tmp)
        for j in 1:nele
            p_x2 = P.P_x2[j]
            ψ_e = tmp[2, j] * p_x2
            ψ_g = tmp[1, j] * p_x2

            T12 = P.P_σ12[j, i - 1]
            T11 = P.P_σ11[j, i - 1]
            T22 = P.P_σ22[j, i - 1]
            T21 = conj(T12)

            ψs[2, j, i] = T11 * ψ_e + T12 * ψ_g
            ψs[1, j, i] = T22 * ψ_g + T21 * ψ_e
        end
    end
    ψs
end
