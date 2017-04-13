#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

# Propagate the state `ϕ` by `nstep` number of time steps each of length `δt`
# The decay matrix is `Γ`, the element `(i,j)` of this matrix represents the decay rate
# from state `j` to state `i`.
# The Rabi flopping happens on resonant from state `i1` to state `i2` with Rabi frequency `Ω`
function propagate(ϕ, δt, nstep, Γ, i1, i2, Ω, nstate, rates, buff1, rd)
    # Do propagation of coherent and dissipative part separately.
    # Use precise exponentiation for the (potentially fast) Rabi flopping and use linear
    # approximation for the decay probability.
    # Compute the decay probability in the middle of the time step to improve precision.
    if nstep == 0
        return ϕ
    end

    # Half step
    @fastmath sinΩt_2 = sin(Ω * δt / 2)
    @fastmath cosΩt_2 = cos(Ω * δt / 2)
    # Full step
    @fastmath sinΩt = sin(Ω * δt)
    @fastmath cosΩt = cos(Ω * δt)
    # Decay rate
    T = real(eltype(ϕ))
    s::T = 0
    for i in 1:nstate
        s += abs2(ϕ[i])
    end

    # First compute the shifted wavefunction. We'll propagate full steps from here
    # until the last step.
    # This also checks bound so we can skip bound checking later
    ϕ1 = ϕ[i1]
    ϕ2 = ϕ[i2]
    ϕ1, ϕ2 = ϕ1 * cosΩt_2 + ϕ2 * sinΩt_2, ϕ2 * cosΩt_2 - ϕ1 * sinΩt_2
    ϕ[i1] = ϕ1
    ϕ[i2] = ϕ2
    @inbounds for i in 1:(nstep - 1)
        # We are at `t = (i - 1/2) * δt`
        # Decide if we need to decay first
        @simd for j in 1:nstate
            buff1[j] = abs2(ϕ[j]) * rates[j] * δt
        end
        # buff1 .= abs2.(ϕ) .* rates .* δt
        r = rand(rd) * s
        decay = 0
        for j in 1:nstate
            r -= buff1[j]
            if r < 0
                decay = j
                break
            end
        end
        if decay != 0
            # decay happend from state `decay`
            # Decide which state it'll fall into
            r = rand(rd) * rates[decay]
            local k = 0
            for k in 1:nstate
                r -= Γ[k, decay]
                if r < 0
                    break
                end
            end
            s = 1
            ϕ .= 0
            if i1 == k
                ϕ[i1] = cosΩt_2
                ϕ[i2] = -sinΩt_2
            elseif i2 == k
                ϕ[i1] = sinΩt_2
                ϕ[i2] = cosΩt_2
            else
                ϕ[k] = 1
            end
            continue
        end
        # No decay, cool. Propagate by a full time step to `t = (i + 1/2) * δt`
        # First propagate the decay term
        @fastmath factor = 1 / 2 / sqrt(s)
        s = 0 # normalization
        @simd for j in 1:nstate
            v = ϕ[j] * (1 - rates[j] * factor * δt)
            s += abs2(v)
            ϕ[j] = v
        end
        # Then propagate the drive
        ϕ1 = ϕ[i1]
        ϕ2 = ϕ[i2]
        ϕ1, ϕ2 = ϕ1 * cosΩt + ϕ2 * sinΩt, ϕ2 * cosΩt - ϕ1 * sinΩt
        ϕ[i1] = ϕ1
        ϕ[i2] = ϕ2
    end
    ϕ1 = ϕ[i1]
    ϕ2 = ϕ[i2]
    ϕ1, ϕ2 = ϕ1 * cosΩt_2 + ϕ2 * sinΩt_2, ϕ2 * cosΩt_2 - ϕ1 * sinΩt_2
    ϕ[i1] = ϕ1
    ϕ[i2] = ϕ2
    ϕ ./= √(s)
    return ϕ
end

function average(ϕ0, δt, nstep, Γ, i1, i2, Ω, n, rates, rd)
    len = length(ϕ0)
    ϕ = similar(ϕ0)
    sumϕ = zeros(real(eltype(ϕ0)), len)
    sumϕ² = zeros(real(eltype(ϕ0)), len)
    T = real(eltype(ϕ))
    buff1 = Vector{T}(len)
    @inbounds for i in 1:n
        ϕ .= ϕ0
        propagate(ϕ, δt, nstep, Γ, i1, i2, Ω, len, rates, buff1, rd)
        for j in 1:len
            ϕi = abs2(ϕ[j])
            sumϕ[j] += ϕi
            sumϕ²[j] += ϕi^2
        end
    end
    @inbounds for i in 1:len
        avgϕi = sumϕ[i] / n
        avgϕ²i = sumϕ²[i] / n
        sumϕ[i] = avgϕi
        sumϕ²[i] = (avgϕ²i - avgϕi^2) / sqrt(n - 1)
    end
    return sumϕ, sumϕ²
end


ϕ = [1.0, 0.0]
δt = 1e-7
nstep = 100_000
Γ = [0 3e4
      2e4 0]
i1 = 1
i2 = 2
Ω = 2π * 400e3

using PyPlot
pts = 0:10:1000
res = Vector{Float64}(length(pts))
unc = Vector{Float64}(length(pts))

function f(pts, ϕ, Γ, i1, i2, Ω, res, unc)
    len = length(ϕ)
    rates = Vector{eltype(Γ)}(len)
    for i in 1:len
        local s = zero(eltype(Γ))
        for j in 1:len
            s += Γ[j, i]
        end
        rates[i] = s
    end
    rds = [MersenneTwister(0) for i in 1:Threads.nthreads()]
    # @show average(ϕ, 1.1e-7, 0, Γ, i1, i2, Ω, 10000, rates, rds[Threads.threadid()])
    @time Threads.@threads for i in 1:length(pts)
        local a, s
        a, s = average(ϕ, 1.1e-7, pts[i], Γ, i1, i2, Ω, 10000, rates, rds[Threads.threadid()])
        res[i] = a[1]
        unc[i] = s[1]
    end
end
f(pts, ϕ, Γ, i1, i2, Ω, res, unc)
errorbar(pts, res, unc)
# # for i in 1:npts
# #     a, s = average(ϕ, 1.1e-7, 10 * (i - 1), Γ, i1, i2, 0, 10000)
# #     res[i] = a[1]
# #     unc[i] = s[1]
# # end
# # errorbar(1:npts, res, unc)

grid()
show()

# @time @show average(ϕ, 5.5e-7, 200, Γ, i1, i2, 0, 1)
# @time @show average(ϕ, 2.75e-7, 40_000, Γ, i1, i2, Ω, 100)
# @time @show average(ϕ, 1.1e-7, 100_000, Γ, i1, i2, Ω, 100)
