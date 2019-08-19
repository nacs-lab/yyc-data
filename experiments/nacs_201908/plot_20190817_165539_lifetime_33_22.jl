#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190817_165539.mat",
                "data_20190817_233758.mat",
                "data_20190803_125057.mat",
                "data_20190803_002700.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 ]
const specs = [([0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100]),
               ([0.0, 20, 50, 100, 200, 500],
                [0.0, 20, 50, 100, 200, 500],
                [0.0, 10, 20, 50, 100, 200, 500],
                [0.0, 5, 10, 20, 50, 100, 200],
                [0.0, 1, 2, 5, 10, 20, 50, 100]),
               ([0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100],
                [0.0, 1, 2, 5, 10, 20, 50, 100]),
               ([1.0, 2, 3, 4, 5, 8, 10, 20, 50],
                [1.0, 2, 3, 4, 5, 8, 10, 20, 50])]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

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

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

model_exp_off(x, p) = p[1] .* exp.(x .* -p[2]) .+ p[3]
model_exp(x, p) = p[1] .* exp.(x .* -p[2])

data_nacs = [[datas_nacs[1][i]; datas_nacs[2][i]] for i in 1:length(datas_nacs[1])]

data_nacs[5] = [data_nacs[5]; datas_nacs[3][1]; datas_nacs[3][2];
                datas_nacs[3][3]; datas_nacs[3][4]; datas_nacs[3][5];
                datas_nacs[4][1]]

fits_nacs = [(i > 2 ? fit_survival(model_exp_off, data_nacs[i], [0.25, 1 / 20, 0.05]) :
              fit_survival(model_exp, data_nacs[i], [0.25, 1 / 20]))
             for i in 1:length(data_nacs)]

const prefix = joinpath(@__DIR__, "imgs", "data_20190817_165539_lifetime_33_22")

tps = [1, 2, 6, 9, 15];
rs = Float64[]
rs_s = Float64[]

figure()
for i in 1:5
    NaCsPlot.plot_survival_data(data_nacs[i], fmt="C$(i - 1).", label="$(tps[i]) mW")
    plot(fits_nacs[i].plotx, fits_nacs[i].ploty, "C$(i - 1)-")
    push!(rs, fits_nacs[i].param[2] * 1000)
    push!(rs_s, fits_nacs[i].unc[2] * 1000)
end
legend(ncol=2, fontsize="small")
grid()
xlim([0, 510])
title("3, 3 + 2, 2 Lifetime")
xlabel("Hold time (ms)")
ylabel("Two body survival")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

model_pow_off(x, p) = p[1] .* x.^p[2] .+ p[3]
function gen_pow_off_model(pow)
    model(x, p) = p[1] .* x.^pow .+ p[2]
end
fit_rs = fit_data(model_pow_off, tps, rs, rs_s, [0.07, 2.5, 1.0], plot_lo=0)
fit_rs_pm = fit_data(gen_pow_off_model(1.8), tps, rs, rs_s, [0.1, 1.0], plot_lo=0)
fit_rs_pp = fit_data(gen_pow_off_model(3.2), tps, rs, rs_s, [0.1, 1.0], plot_lo=0)

figure()
errorbar(tps, rs, rs_s, fmt="C0.")
plot(fit_rs.plotx, fit_rs.ploty, "C0")
plot(fit_rs_pm.plotx, fit_rs_pm.ploty, "C1--")
plot(fit_rs_pp.plotx, fit_rs_pp.ploty, "C1--")
text(0.3, 50, "\$r=r_0+a\\cdot P^b\$", color="C0")
text(2.2, 23, ("\$r_0=$(fit_rs.uncs[3]) s^{-1}\$\n" *
               "\$a=$(fit_rs.uncs[1])\$\n" *
               "\$b=$(fit_rs.uncs[2])\$"), color="C0")
grid()
xlim([0, 16])
ylim([0, 70])
title("3, 3 + 2, 2 loss rates")
xlabel("Tweezer Power (mW)")
ylabel("Loss rate (\$s^{-1}\$)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_rates")

NaCsPlot.maybe_show()
