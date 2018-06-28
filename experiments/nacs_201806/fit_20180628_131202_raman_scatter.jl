#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Atomic: all_scatter_D

const δf1 = -76.82e9
const δf2 = -76.82e9 - 1.77e9
const rlof_f1 = (61.542e6 / (δf1 - 1.107266e9))^2
const rlof_f2 = (61.542e6 / (δf2 - 1.107266e9))^2
const rhif_f1 = (61.542e6 / (δf1 + 664.360e6))^2
const rhif_f2 = (61.542e6 / (δf2 + 664.360e6))^2

const rates_f1_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f1_down = rates_f1_up
const rates_f1_diagonal = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f1, rlof_f1)
const rates_f2_diagonal = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f2, rlof_f2)
const rates_f2_counterop = all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_f2, rlof_f2)

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

function propagate_f1(A, init, t)
    res = exp(A .* t) * init # expm!!!
    return res[6] + res[7] + res[8]
end

function gen_model(rates, init)
    A = rates_to_A(rates)
    t->propagate_f1(A, init, t)
end

const iname_a = joinpath(@__DIR__, "data", "data_20180628_131202.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :f1_diagonal=>[1, 2, 4, 8, 16, 32, 64, 128, 256],
    :f2_diagonal=>[1, 2, 4, 8, 16, 32, 64, 128, 256],
    :f1_up=>[1, 2, 4, 8, 16, 32, 64, 128, 256],
    :f2_counterop=>[1, 2, 4, 8, 16, 32, 64, 128, 256],
    :f1_down=>[1, 2, 4, 8, 16, 32, 64, 128, 256],
    :t0=>[0],
)

const split_a = NaCsData.split_data(data_a, spec)

const data_f1_diagonal = [split_a[:t0]; split_a[:f1_diagonal]]
const data_f2_diagonal = [split_a[:t0]; split_a[:f2_diagonal]]
const data_f1_up = [split_a[:t0]; split_a[:f1_up]]
const data_f2_counterop = [split_a[:t0]; split_a[:f2_counterop]]
const data_f1_down = [split_a[:t0]; split_a[:f1_down]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180628_131202_raman_scatter")

function gen_fit_model(rates, t_prescale)
    f = gen_model(rates, Float64[1, 0, 0, 0, 0, 0, 0, 0])
    function (x, p)
        p[1] .* f.(x .* p[2] .* t_prescale)
    end
end

function fit_and_plot_op(rates, data)
    ts, _ratios, _uncs = NaCsData.get_values(data)
    ratios = _ratios[:, 2]
    uncs = _uncs[:, 2]
    init = Float64[1, 0, 0, 0, 0, 0, 0, 0]
    init_f1 = Float64[0, 0, 0, 0, 0, 1, 0, 0]
    r = sum(rates * init)
    r_f1 = sum(rates * init_f1)
    τ = 1 / r
    τ_f1 = 1 / r_f1
    τ_max = max(τ, τ_f1)
    tmax = maximum(ts)
    t_prescale = 10τ_max / tmax
    model = gen_fit_model(rates, t_prescale)
    pinit = [1.0, 1.0]
    fit = curve_fit(model, ts, ratios, pinit)
    plot_ts = linspace(0, tmax, 1000)
    ys = model(plot_ts, fit.param)
    plot(plot_ts, ys, color="C1")
    errorbar(ts, ratios, uncs, fmt="C0o")
    grid()
    ylim(0, ylim()[2])
    xlim(0, xlim()[2])
    xlabel("\$t (ms)\$")
    nrm, tscale = fit.param
    nrm_s, tscale_s = estimate_errors(fit)
    tscale *= t_prescale
    tscale_s *= t_prescale
    τ_f1_real = τ_f1 / tscale
    τ_f1_real_s = tscale_s / tscale * τ_f1_real
    τ_f2_real = τ / tscale
    τ_f2_real_s = tscale_s / tscale * τ_f2_real
    return (nrm=Unc(nrm, nrm_s),
            tscale=Unc(tscale, tscale_s),
            τ_f1=Unc(τ_f1_real, τ_f1_real_s),
            τ_f2=Unc(τ_f2_real, τ_f2_real_s))
end

function fit_text(x, y, fit)
    text(x, y, "\$\\tau_{F1}=$(fit.τ_f1)\$ms\n\$\\tau_{F2}=$(fit.τ_f2)\$ms\n\$Total=$(fit.nrm)\$")
end

figure()
fit_f1_diagonal = fit_and_plot_op(rates_f1_diagonal, data_f1_diagonal)
fit_text(100, 0.07, fit_f1_diagonal)
title("F1 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f1_diagonal")

figure()
fit_f2_diagonal = fit_and_plot_op(rates_f2_diagonal, data_f2_diagonal)
fit_text(100, 0.07, fit_f2_diagonal)
title("F2 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f2_diagonal")

figure()
fit_f1_up = fit_and_plot_op(rates_f1_up, data_f1_up)
fit_text(100, 0.09, fit_f1_up)
title("F1 up")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f1_up")

figure()
fit_f2_counterop = fit_and_plot_op(rates_f2_counterop, data_f2_counterop)
fit_text(30, 0.04, fit_f2_counterop)
title("F2 counter-OP")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f2_counterop")

figure()
fit_f1_down = fit_and_plot_op(rates_f1_down, data_f1_down)
fit_text(100, 0.02, fit_f1_down)
title("F1 down")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f1_down")

NaCsPlot.maybe_show()
