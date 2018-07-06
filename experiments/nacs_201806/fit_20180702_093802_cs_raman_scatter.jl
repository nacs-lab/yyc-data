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

const δf3 = -49.229e9 + 9.192631770e9
const δf4 = -49.229e9
const rlof_f3 = (61.542e6 / (δf3 - 5.170855370625e9))^2
const rlof_f4 = (61.542e6 / (δf4 - 5.170855370625e9))^2
const rhif_f3 = (61.542e6 / (δf3 + 4.021776399375e9))^2
const rhif_f4 = (61.542e6 / (δf4 + 4.021776399375e9))^2

const rates_f4_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f4, rlof_f4)
const rates_f4_down = rates_f4_up
const rates_f3_diagonal = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f3, rlof_f3)
const rates_f4_diagonal = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f4, rlof_f4)
const rates_f3_counterop = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f3, rlof_f3)

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

function propagate_f3(A, init, t)
    res = exp(A .* t) * init # expm!!!
    return res[6] + res[7] + res[8]
end

function gen_model(rates, init)
    A = rates_to_A(rates)
    t->propagate_f3(A, init, t)
end

const iname_a = joinpath(@__DIR__, "data", "data_20180702_093802.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[4, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :f3_diagonal=>[4, 8, 16, 32, 64, 128, 256],
    :f4_diagonal=>[4, 8, 16, 32, 64, 128, 256],
    :f3_counterop=>[4, 8, 16, 32, 64, 128, 256],
    :f4_up=>[4, 8, 16, 32, 64, 128, 256],
    :f4_down=>[4, 8, 16, 32, 64, 128, 256],
    :bg=>[4, 8, 16, 32, 64, 128, 256],
    :t0=>[0],
)

const split_a = NaCsData.split_data(data_a, spec)

const data_f3_diagonal = [split_a[:t0]; split_a[:f3_diagonal]]
const data_f4_diagonal = [split_a[:t0]; split_a[:f4_diagonal]]
const data_f4_up = [split_a[:t0]; split_a[:f4_up]]
const data_f3_counterop = [split_a[:t0]; split_a[:f3_counterop]]
const data_f4_down = [split_a[:t0]; split_a[:f4_down]]
const data_bg = [split_a[:t0]; split_a[:bg]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180702_093802_cs_raman_scatter")

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
    init_f3 = Float64[0, 0, 0, 0, 0, 1, 0, 0]
    r = sum(rates * init)
    r_f3 = sum(rates * init_f3)
    τ = 1 / r
    τ_f3 = 1 / r_f3
    τ_max = max(τ, τ_f3)
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
    τ_f3_real = τ_f3 / tscale
    τ_f3_real_s = tscale_s / tscale * τ_f3_real
    τ_f4_real = τ / tscale
    τ_f4_real_s = tscale_s / tscale * τ_f4_real
    return (nrm=Unc(nrm, nrm_s),
            tscale=Unc(tscale, tscale_s),
            τ_f3=Unc(τ_f3_real, τ_f3_real_s),
            τ_f4=Unc(τ_f4_real, τ_f4_real_s))
end

function fit_text(x, y, fit)
    text(x, y, "\$\\tau_{F3}=$(fit.τ_f3)\$ms\n\$\\tau_{F4}=$(fit.τ_f4)\$ms\n\$Total=$(fit.nrm)\$")
end

figure()
fit_f3_diagonal = fit_and_plot_op(rates_f3_diagonal, data_f3_diagonal)
fit_text(100, 0.07, fit_f3_diagonal)
title("F3 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f3_diagonal")

figure()
fit_f4_diagonal = fit_and_plot_op(rates_f4_diagonal, data_f4_diagonal)
fit_text(100, 0.07, fit_f4_diagonal)
title("F4 diagonal")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f4_diagonal")

figure()
fit_f4_up = fit_and_plot_op(rates_f4_up, data_f4_up)
fit_text(100, 0.09, fit_f4_up)
title("F3 up")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f4_up")

figure()
fit_f3_counterop = fit_and_plot_op(rates_f3_counterop, data_f3_counterop)
fit_text(30, 0.04, fit_f3_counterop)
title("F4 counter-OP")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f3_counterop")

figure()
fit_f4_down = fit_and_plot_op(rates_f4_down, data_f4_down)
fit_text(100, 0.02, fit_f4_down)
title("F3 down")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f4_down")

NaCsPlot.maybe_show()
