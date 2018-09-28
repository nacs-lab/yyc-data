#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Atomic: all_scatter_D, all_scatter_bothD

const δcs_f3 = -45.7e9 + 9.192631770e9
const δcs_f4 = -45.7e9
const δcs_tw_d2 = -63.296e12
const δcs_tw_d1 = -66.686e12
const rlof_cs_f3 = (32.889e6 / (δcs_f3 - 5.170855370625e9))^2
const rlof_cs_f4 = (32.889e6 / (δcs_f4 - 5.170855370625e9))^2
const rhif_cs_f3 = (32.889e6 / (δcs_f3 + 4.021776399375e9))^2
const rhif_cs_f4 = (32.889e6 / (δcs_f4 + 4.021776399375e9))^2

const rlof_cs_tw_d2 = (32.889e6 / (δcs_tw_d2 - 5.170855370625e9))^2
const rlof_cs_tw_d1 = (28.743e6 / (δcs_tw_d1 - 5.170855370625e9))^2
const rhif_cs_tw_d2 = (32.889e6 / (δcs_tw_d2 + 4.021776399375e9))^2
const rhif_cs_tw_d1 = (28.743e6 / (δcs_tw_d1 + 4.021776399375e9))^2

const rates_cs_f4_up = all_scatter_D(true, 7, (0.25, 0.5, 0.25), rhif_cs_f4, rlof_cs_f4)
const rates_cs_f4_down = rates_cs_f4_up
const rates_cs_f3_diagonal = all_scatter_D(true, 7, (0.5, 0.0, 0.5), rhif_cs_f3, rlof_cs_f3)
const rates_cs_f4_diagonal = all_scatter_D(true, 7, (0.25, 0.5, 0.25), rhif_cs_f4, rlof_cs_f4)
const rates_cs_f3_counterop = all_scatter_D(true, 7, (0.5, 0.0, 0.5), rhif_cs_f3, rlof_cs_f3)
const rates_cs_tw = all_scatter_bothD(7, (0.0, 1.0, 0.0)) do D2, Fgx2, mFgx2, Fex2, mFex2
    if D2
        return Fgx2 == 8 ? sqrt(rhif_cs_tw_d2) : sqrt(rlof_cs_tw_d2)
    else
        return Fgx2 == 8 ? sqrt(rhif_cs_tw_d1) : sqrt(rlof_cs_tw_d1)
    end
end

const δna_f1 = -76.82e9
const δna_f2 = -76.82e9 - 1.77e9
const rlof_na_f1 = (61.542e6 / (δna_f1 - 1.107266e9))^2
const rlof_na_f2 = (61.542e6 / (δna_f2 - 1.107266e9))^2
const rhif_na_f1 = (61.542e6 / (δna_f1 + 664.360e6))^2
const rhif_na_f2 = (61.542e6 / (δna_f2 + 664.360e6))^2

const rates_na_f1_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_na_f1, rlof_na_f1)
const rates_na_f1_down = rates_na_f1_up
const rates_na_f1_diagonal = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_na_f1, rlof_na_f1)
const rates_na_f2_diagonal = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_na_f2, rlof_na_f2)
const rates_na_f2_counterop = all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_na_f2, rlof_na_f2)

function rates_to_A(rates)
    nx, ny = size(rates)
    A = Matrix{Float64}(nx, ny)
    @inbounds for i in 1:nx
        s = 0.0
        for j in 1:ny
            r = rates[j, i]
            A[j, i] = r
            s += r
        end
        A[i, i] -= s
    end
    return A
end

function propagate_lof(A, init, t)
    res = exp(A .* t) * init # expm!!!
    n_2 = length(res) ÷ 2
    return sum(res[n_2 + 2:end])
end

function gen_model(rates, init)
    A = rates_to_A(rates)
    t->propagate_lof(A, init, t)
end

function gen_model2(rates1, rates2, init)
    A1 = rates_to_A(rates1)
    A2 = rates_to_A(rates2)
    (t, a2)->propagate_lof(A1 .+ a2 .* A2, init, t)
end

function gen_fit_model(rates, t_prescale, τ, p0)
    init = zeros(Float64, size(rates, 1))
    init[1] = 1
    f = gen_model(rates, init)
    function (x, p)
        p0 .* f.(x .* p[1] .* t_prescale) .* exp.(.- x ./ τ)
    end
end

function gen_fit_model′(rates, t_prescale, τ)
    init = zeros(Float64, size(rates, 1))
    init[1] = 1
    f = gen_model(rates, init)
    function (x, p)
        p[2] .* f.(x .* p[1] .* t_prescale) .* exp.(.- x ./ τ)
    end
end

function gen_fit_model2(rates1, rates2, t_prescale, τ, p0)
    init = zeros(Float64, size(rates1, 1))
    init[1] = 1
    f = gen_model2(rates1, rates2, init)
    function (x, p)
        p0 .* f.(x .* t_prescale, p[1]) .* exp.(.- x ./ τ)
    end
end

function gen_fit_model2′(rates1, rates2, t_prescale, τ)
    init = zeros(Float64, size(rates1, 1))
    init[1] = 1
    f = gen_model2(rates1, rates2, init)
    function (x, p)
        p[2] .* f.(x .* t_prescale, p[1]) .* exp.(.- x ./ τ)
    end
end

function fit_and_plot_op(rates, data, τlifetime; p0=nothing, rates0=nothing, kws...)
    ts, _ratios, _uncs = NaCsData.get_values(data)
    ratios = _ratios[:, 2]
    uncs = _uncs[:, 2]
    init = zeros(Float64, size(rates, 1))
    init[1] = 1
    init_lof = zeros(Float64, size(rates, 1))
    init_lof[size(rates, 1) ÷ 2 + 2] = 1
    r = sum(rates * init)
    r_lof = sum(rates * init_lof)
    τ = 1 / r
    τ_lof = 1 / r_lof
    τ_max = max(τ, τ_lof)
    tmax = maximum(ts)
    t_prescale = 10τ_max / tmax
    if p0 === nothing
        model = (rates0 === nothing ? gen_fit_model′(rates, t_prescale, τlifetime) :
                 gen_fit_model2′(rates0, rates, t_prescale, τlifetime))
        pinit = [1.0, 1.0]
    else
        model = (rates0 === nothing ? gen_fit_model(rates, t_prescale, τlifetime, p0.a) :
                 gen_fit_model2(rates0, rates, t_prescale, τlifetime, p0.a))
        pinit = [1.0]
    end
    fit = curve_fit(model, ts, ratios, pinit)
    plot_ts = linspace(0, tmax, 1000)
    ys = model(plot_ts, fit.param)
    plot(plot_ts, ys, color="C0")
    errorbar(ts, ratios, uncs; fmt="C0.", kws...)
    grid()
    ylim(0, ylim()[2])
    xlim(0, xlim()[2])
    xlabel("\$t (ms)\$")
    if p0 === nothing
        tscale, nrm = fit.param
        tscale_s, nrm_s = estimate_errors(fit)
        nrm_u = Unc(nrm, nrm_s)
    else
        tscale, = fit.param
        tscale_s, = estimate_errors(fit)
        nrm_u = p0
    end
    tscale *= t_prescale
    tscale_s *= t_prescale
    τ_lof_real = τ_lof / tscale
    τ_lof_real_s = tscale_s / tscale * τ_lof_real
    τ_hif_real = τ / tscale
    τ_hif_real_s = tscale_s / tscale * τ_hif_real
    return (nrm=nrm_u,
            tscale=Unc(tscale, tscale_s),
            τ_lof=Unc(τ_lof_real, τ_lof_real_s),
            τ_hif=Unc(τ_hif_real, τ_hif_real_s))
end

function fit_survival(model, data, p0; plotx=nothing, use_unc=false, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        plotx = linspace(mid - span * plot_scale / 2, mid + span * plot_scale / 2, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2].^2, p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

lifetime_model(x, p) = p[1] .* exp.(.-x ./ p[2])

const iname_a = joinpath(@__DIR__, "data", "data_20180903_213034.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec = OrderedDict(
    :lo_diagonal=>[1, 2, 4, 8, 16, 32, 64, 128, 256, 512],
    :hi_diagonal=>[1, 2, 4, 8, 16, 32, 64, 128, 256, 512],
    :up=>[1, 2, 4, 8, 16, 32, 64, 128, 256, 512],
    :counterop=>[1, 2, 4, 8, 16, 32, 64, 128, 256, 512],
    :down=>[1, 2, 4, 8, 16, 32, 64, 128, 256, 512],
    :t0=>[0],
    :tw=>[16, 32, 64, 128, 256, 512, 1024, 2048],
    :lifetime=>[16, 32, 64, 128, 256, 512, 1024, 2048],
)

const split_na_a = NaCsData.split_data(data_na_a, spec)
const split_cs_a = NaCsData.split_data(data_cs_a, spec)

const data_cs_lo_diagonal = [split_cs_a[:t0]; split_cs_a[:lo_diagonal]]
const data_cs_hi_diagonal = [split_cs_a[:t0]; split_cs_a[:hi_diagonal]]
const data_cs_up = [split_cs_a[:t0]; split_cs_a[:up]]
const data_cs_counterop = [split_cs_a[:t0]; split_cs_a[:counterop]]
const data_cs_down = [split_cs_a[:t0]; split_cs_a[:down]]
const data_cs_tw = [split_cs_a[:t0]; split_cs_a[:tw]]
const data_cs_lifetime = split_cs_a[:lifetime]

const data_na_lo_diagonal = [split_na_a[:t0]; split_na_a[:lo_diagonal]]
const data_na_hi_diagonal = [split_na_a[:t0]; split_na_a[:hi_diagonal]]
const data_na_up = [split_na_a[:t0]; split_na_a[:up]]
const data_na_counterop = [split_na_a[:t0]; split_na_a[:counterop]]
const data_na_down = [split_na_a[:t0]; split_na_a[:down]]
const data_na_tw = [split_na_a[:t0]; split_na_a[:tw]]
const data_na_lifetime = split_na_a[:lifetime]

const prefix = joinpath(@__DIR__, "imgs", "data_20180903_213034_raman_scatter")

fit_cs = fit_survival(lifetime_model, data_cs_lifetime, [0.95, 4.4],
                      plotx=linspace(0, 2500, 1001))
fit_na = fit_survival(lifetime_model, data_na_lifetime, [0.95, 4.4],
                      plotx=linspace(0, 2500, 1001))

figure()
NaCsPlot.plot_survival_data(data_cs_lifetime, fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-", label="Cs")
NaCsPlot.plot_survival_data(data_na_lifetime, fmt="C1.")
plot(fit_na.plotx, fit_na.ploty, "C1-", label="Na")
text(100, 0.55, "\$\\tau_{Na}=$(Unc(fit_na.param[2] / 1000, fit_na.unc[2] / 1000))s\$", color="C0")
text(100, 0.45, "\$\\tau_{Cs}=$(Unc(fit_cs.param[2] / 1000, fit_cs.unc[2] / 1000))s\$", color="C1")
grid()
legend()
xlim([0, 2500])
ylim([0.4, 1])
title("Lifetime")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_lifetime")

const p0_na = Unc(fit_na.param[1], fit_na.unc[1])
const p0_cs = Unc(fit_cs.param[1], fit_cs.unc[1])
const τ_na = fit_na.param[2]
const τ_cs = fit_cs.param[2]

function fit_text_cs(x, y, fit; kws...)
    text(x, y, "\$\\tau_{F3}=$(fit.τ_lof)\$ms\n\$\\tau_{F4}=$(fit.τ_hif)\$ms\n\$Total=$(fit.nrm)\$"; kws...)
end

function fit_text_na(x, y, fit; kws...)
    text(x, y, "\$\\tau_{F1}=$(fit.τ_lof)\$ms\n\$\\tau_{F2}=$(fit.τ_hif)\$ms\n\$Total=$(fit.nrm)\$"; kws...)
end

figure()
fit_cs_tw = fit_and_plot_op(rates_cs_tw, data_cs_tw, τ_cs, p0=p0_cs, label="Cs")
NaCsPlot.plot_survival_data(data_na_tw, fmt="C1.-", label="Na")
fit_text_cs(1000, 0.06, fit_cs_tw, color="C0")
legend()
title("Tweezer")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_tw")
rates_cs0 = rates_cs_tw ./ fit_cs_tw.tscale.a

figure()
fit_cs_f3_diagonal = fit_and_plot_op(rates_cs_f3_diagonal, data_cs_lo_diagonal,
                                     τ_cs, rates0=rates_cs0)
fit_text_cs(100, 0.07, fit_cs_f3_diagonal)
title("Cs F3 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_f3_diagonal")

figure()
fit_cs_f4_diagonal = fit_and_plot_op(rates_cs_f4_diagonal, data_cs_hi_diagonal,
                                     τ_cs, rates0=rates_cs0)
fit_text_cs(100, 0.07, fit_cs_f4_diagonal)
title("Cs F4 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_f4_diagonal")

figure()
fit_cs_f4_up = fit_and_plot_op(rates_cs_f4_up, data_cs_up,
                               τ_cs, rates0=rates_cs0)
fit_text_cs(100, 0.09, fit_cs_f4_up)
title("Cs F4 up")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_f4_up")

figure()
fit_cs_f3_counterop = fit_and_plot_op(rates_cs_f3_counterop, data_cs_counterop,
                                      τ_cs, rates0=rates_cs0)
fit_text_cs(30, 0.04, fit_cs_f3_counterop)
title("Cs F3 counter-OP")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_f3_counterop")

figure()
fit_cs_f4_down = fit_and_plot_op(rates_cs_f4_down, data_cs_down,
                                 τ_cs, rates0=rates_cs0)
fit_text_cs(100, 0.02, fit_cs_f4_down)
title("Cs F4 down")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_f4_down")

figure()
fit_na_f1_diagonal = fit_and_plot_op(rates_na_f1_diagonal, data_na_lo_diagonal, τ_na, p0=p0_na)
fit_text_na(100, 0.07, fit_na_f1_diagonal)
title("Na F1 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_f1_diagonal")

figure()
fit_na_f2_diagonal = fit_and_plot_op(rates_na_f2_diagonal, data_na_hi_diagonal, τ_na, p0=p0_na)
fit_text_na(100, 0.07, fit_na_f2_diagonal)
title("Na F2 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_f2_diagonal")

figure()
fit_na_f1_up = fit_and_plot_op(rates_na_f1_up, data_na_up, τ_na, p0=p0_na)
fit_text_na(100, 0.09, fit_na_f1_up)
title("Na F1 up")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_f1_up")

figure()
fit_na_f2_counterop = fit_and_plot_op(rates_na_f2_counterop, data_na_counterop, τ_na, p0=p0_na)
fit_text_na(15, 0.01, fit_na_f2_counterop)
title("Na F2 counter-OP")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_f2_counterop")

figure()
fit_na_f1_down = fit_and_plot_op(rates_na_f1_down, data_na_down, τ_na, p0=p0_na)
fit_text_na(100, 0.02, fit_na_f1_down)
title("Na F1 down")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_f1_down")

NaCsPlot.maybe_show()