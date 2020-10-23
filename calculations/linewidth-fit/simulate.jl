#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive, binomial_estimate
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using Distributions
using Statistics

function generate_data(f, f0, fwhm, pmin, pmax, trail)
    p0 = 1 - 1 / (1 + (2 * (f - f0) / fwhm)^2)
    p = pmin + (pmax - pmin) * p0
    return binomial_estimate(rand(Binomial(trail, p)), trail)
end

function generate_data(fs::AbstractArray, f0, fwhm, pmin, pmax, trail)
    nf = length(fs)
    res = Vector{Float64}(undef, nf)
    unc = Vector{Float64}(undef, nf)
    for i in 1:nf
        res[i], unc[i] = generate_data(fs[i], f0, fwhm, pmin, pmax, trail)
    end
    return res, unc
end

function model(x, p)
    f0, fwhm, pmin, pmax = p
    return pmin .+ (pmax - pmin) .* (1 .- 1 ./ (1 .+ (2 .* (x .- f0) ./ fwhm).^2))
end

function generate_fit(fs::AbstractArray, f0, fwhm, pmin, pmax, trail)
    data = generate_data(fs, f0, fwhm, pmin, pmax, trail)
    return fit_data(model, fs, data[1], data[2], [f0, fwhm, pmin, pmax])
end

function generate_fits(fs::AbstractArray, f0, fwhm, pmin, pmax, trail, trail_fit)
    vals = Vector{Float64}(undef, trail_fit)
    uncs = Vector{Float64}(undef, trail_fit)
    for i in 1:trail_fit
        @label retry
        fit = try
            generate_fit(fs, f0, fwhm, pmin, pmax, trail)
        catch
            @goto retry
        end
        width = fit.param[2]
        width_unc = fit.unc[2]
        if width < 0 || width_unc >= width || width < 0.1 * fwhm || width > 10 * fwhm
            @goto retry
        end
        vals[i] = width
        uncs[i] = width_unc
    end
    return mean(vals), std(vals), mean(uncs), std(uncs)
end

function generate_fits_trails(fs::AbstractArray, f0, fwhm, pmin, pmax,
                              trails::AbstractArray, trail_fit)
    nt = length(trails)
    vs = Vector{Float64}(undef, nt)
    v_uncs = Vector{Float64}(undef, nt)
    us = Vector{Float64}(undef, nt)
    u_uncs = Vector{Float64}(undef, nt)
    for i in 1:nt
        @show trails[i]
        vs[i], v_uncs[i], us[i], u_uncs[i] = generate_fits(fs, f0, fwhm, pmin, pmax,
                                                           trails[i], trail_fit)
    end
    return vs, v_uncs, us, u_uncs
end

const prefix = joinpath(@__DIR__, "imgs", "sim")

trails = [50, 100, 200, 300, 500, 700, 1000]

freqs = [-1000, -60, -45, -30, -15, 0, 15, 30, 1000] .+ 7
res = generate_fits_trails(freqs, -15, 30, 0.5, 0.8, trails, 1000)
figure(figsize=[12.8, 4.8])
subplot(121)
errorbar(trails, res[1], res[2], label="Value")
errorbar(trails .+ 2, res[3], res[4], label="Uncertainty")
legend()
title("Fitting Result")
xlabel("Trails")
grid()

subplot(122)
plot(trails, (res[1] .- 30) ./ res[3])
title("Relative Bias")
xlabel("Trails")
grid()
NaCsPlot.maybe_save("$(prefix)_good15_7")

freqs = [-1000, -60, -45, -30, -15, 0, 15, 30, 1000]
res = generate_fits_trails(freqs, -15, 30, 0.5, 0.8, trails, 1000)
figure(figsize=[12.8, 4.8])
subplot(121)
errorbar(trails, res[1], res[2], label="Value")
errorbar(trails .+ 2, res[3], res[4], label="Uncertainty")
legend()
title("Fitting Result")
xlabel("Trails")
grid()

subplot(122)
plot(trails, (res[1] .- 30) ./ res[3])
title("Relative Bias")
xlabel("Trails")
grid()
NaCsPlot.maybe_save("$(prefix)_good15")

freqs = [-100, -80, -60, -40, -20, 0.0, 20, 40, 100]
res = generate_fits_trails(freqs, -15, 30, 0.5, 0.8, trails, 1000)
figure(figsize=[12.8, 4.8])
subplot(121)
errorbar(trails, res[1], res[2], label="Value")
errorbar(trails .+ 2, res[3], res[4], label="Uncertainty")
legend()
title("Fitting Result")
xlabel("Trails")
grid()

subplot(122)
plot(trails, (res[1] .- 30) ./ res[3])
title("Relative Bias")
xlabel("Trails")
grid()
NaCsPlot.maybe_save("$(prefix)_bad")

NaCsPlot.maybe_show()
