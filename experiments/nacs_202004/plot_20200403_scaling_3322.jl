#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

function merge_unc(avg, unc)
    D = sum(unc.^-2)
    N = sum(avg .* unc.^-2)
    return N / D, 1 / sqrt(D)
end

function Base.mean(uncs::AbstractArray{T}) where T<:Unc
    return Unc(merge_unc([u.a for u in uncs],
                         [u.s for u in uncs])...)
end

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, params, ratios, uncs, p0;
                  plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
    use_unc = uncs !== nothing
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        if plot_lo === nothing
            plot_lo = mid - span * plot_scale / 2
            if plot_lo * lo <= 0
                plot_lo = 0
            end
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
            if plot_hi * hi <= 0
                plot_hi = 0
            end
        end
        plotx = linspace(plot_lo, plot_hi, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

const powers = [15, 6, 3]
const omega_r = [4.00, 1.178, 0.483] # kHz
const omega_r_unc = [0.10, 0.061, 0.026]
const gamma_m = [1.78, 0.580, 0.280] # kHz
const gamma_m_unc = [0.20, 0.091, 0.048]
const gamma_a = [22.6, 2.37, 0.39] # Hz
const gamma_a_unc = [0.33, 0.50, 0.13]

const omega_m = Unc(91.2, 4.3) # MHz/mW^-1/2
const detuning = 145
const gamma_e = Unc.(gamma_m, gamma_m_unc) .* (4 * detuning^2) ./ (omega_m.^2 .* powers) # GHz
const gamma_e_avg = mean(Unc.(gamma_m, gamma_m_unc) ./ powers) * (4 * detuning^2) / omega_m^2

const omega_a = 4 .* Unc.(omega_r, omega_r_unc) .* detuning ./ omega_m ./ sqrt.(powers) # MHz
const gamma_ea = Unc.(gamma_a, gamma_a_unc) .* (4 * detuning^2) ./ omega_a.^2 ./ 1000 # GHz

const prefix = joinpath(@__DIR__, "imgs", "data_20200403_scaling_3322")

function gen_power_model(pwr)
    return function (x, p)
        return p[1] .* x.^pwr
    end
end
function power_model(x, p)
    return p[1] .* x.^p[2]
end

fit_omega_r = fit_data(gen_power_model(1.375), powers, omega_r, omega_r_unc, [1.0]; plot_lo=2.7)
fit_omega_r_p = fit_data(power_model, powers, omega_r, omega_r_unc, [1.0, 1.375]; plot_lo=2.7)
fit_gamma_m = fit_data(gen_power_model(1), powers, gamma_m, gamma_m_unc, [1.0]; plot_lo=2.7)
fit_gamma_m2 = fit_data(gen_power_model(2), powers, gamma_m, gamma_m_unc, [1.0]; plot_lo=2.7)
fit_gamma_m_p = fit_data(power_model, powers, gamma_m, gamma_m_unc, [1.0, 1]; plot_lo=2.7)
fit_gamma_a = fit_data(gen_power_model(1.75), powers, gamma_a, gamma_a_unc, [1.0]; plot_lo=2.7)
fit_gamma_a2 = fit_data(gen_power_model(2.75), powers, gamma_a, gamma_a_unc, [1.0]; plot_lo=2.7)
fit_gamma_a_p = fit_data(power_model, powers, gamma_a, gamma_a_unc, [1.0, 2.75]; plot_lo=2.7)
fit_omega_a = fit_data(gen_power_model(0.875), powers,
                       [v.a for v in omega_a], [v.s for v in omega_a], [1.0]; plot_lo=2.7)
# Don't fit with power_model on omega_a directly. Correlated uncertainty!
fit_gamma_ea = fit_data(gen_power_model(1), powers,
                        [v.a for v in gamma_ea], [v.s for v in gamma_ea], [1.0]; plot_lo=2.7)
# Don't fit with power_model on gamma_ea directly. Correlated uncertainty!

figure(figsize=[12.6, 16.8])

subplot(3, 2, 1)
errorbar(powers, omega_r, omega_r_unc, fmt="C0s")
plot(fit_omega_r.plotx, fit_omega_r.ploty, "C0", label="\$p^{1.375}\$")
plot(fit_omega_r_p.plotx, fit_omega_r_p.ploty, "C2", label="\$p^n\$")
text(4.5, 0.58, "\$n=$(fit_omega_r_p.uncs[2])\$", color="C2")
legend()
grid()
xscale("log")
yscale("log")
xlabel("Power (mW)")
ylabel("\$\\Omega_{Raman} (2\\pi\\cdot kHz\$)")
title("Raman Rabi frequency")

subplot(3, 2, 2)
errorbar(powers, gamma_m, gamma_m_unc, fmt="C0s")
plot(fit_gamma_m.plotx, fit_gamma_m.ploty, "C0", label="\$p^{1}\$")
plot(fit_gamma_m2.plotx, fit_gamma_m2.ploty, "C1", label="\$p^{2}\$")
plot(fit_gamma_m_p.plotx, fit_gamma_m_p.ploty, "C2", label="\$p^n\$")
text(4.5, 0.30, "\$n=$(fit_gamma_m_p.uncs[2])\$", color="C2")
legend()
grid()
xscale("log")
yscale("log")
xlabel("Power (mW)")
ylabel("\$\\Gamma_{molecule} (2\\pi\\cdot kHz\$)")
title("Molecular scattering rate")

subplot(3, 2, 3)
errorbar(powers, gamma_a, gamma_a_unc, fmt="C0s")
plot(fit_gamma_a.plotx, fit_gamma_a.ploty, "C0", label="\$p^{1.75}\$")
plot(fit_gamma_a2.plotx, fit_gamma_a2.ploty, "C1", label="\$p^{2.75}\$")
plot(fit_gamma_a_p.plotx, fit_gamma_a_p.ploty, "C2", label="\$p^n\$")
text(4.5, 0.30, "\$n=$(fit_gamma_a_p.uncs[2])\$", color="C2")
legend()
grid()
xscale("log")
yscale("log")
xlabel("Power (mW)")
ylabel("\$\\Gamma_{atomic} (2\\pi\\cdot Hz\$)")
title("Atomic scattering (PA) rate")

subplot(3, 2, 4)
errorbar(powers, [v.a for v in gamma_e], [v.s for v in gamma_e], fmt="C0s")
plot([2.7, 15.5], [gamma_e_avg.a, gamma_e_avg.a], "C0")
text(2.55, 1.2, "\$\\Gamma_{excited}=2\\pi\\times$(gamma_e_avg)\\ GHz\$", color="C0")
grid()
xlabel("Power (mW)")
ylabel("\$\\Gamma_{excited} (2\\pi\\cdot GHz\$)")
title("Excited state linewidth\$_{(\\mathrm{from}\\ \\Gamma_{molecule})}\$")

subplot(3, 2, 5)
errorbar(powers, [v.a for v in omega_a], [v.s for v in omega_a], fmt="C0s")
plot(fit_omega_a.plotx, fit_omega_a.ploty, "C0", label="\$p^{0.875}\$")
legend()
grid()
xscale("log")
yscale("log")
xlabel("Power (mW)")
ylabel("\$\\Omega_{atomic} (2\\pi\\cdot MHz\$)")
title("Atomic Rabi frequency \$_{\\mathrm{(both\\ legs)}}\$")

subplot(3, 2, 6)
errorbar(powers, [v.a for v in gamma_ea], [v.s for v in gamma_ea], fmt="C0s")
plot(fit_gamma_ea.plotx, fit_gamma_ea.ploty, "C0", label="\$p^{1}\$")
legend()
grid()
xlabel("Power (mW)")
ylabel("\$\\Gamma_{excited} (2\\pi\\cdot GHz\$)")
title("Excited state linewidth\$_{(\\mathrm{from}\\ \\Gamma_{atomic})}\$")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
