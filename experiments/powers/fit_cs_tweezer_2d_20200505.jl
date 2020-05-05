#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci
using DelimitedFiles
using Printf

const iname_pa2d = joinpath(@__DIR__, "data", "cs_tweezer_pa2d_20200505.csv")
const data_pa2d = readdlm(iname_pa2d, ',', Float64, skipstart=1)

data_pa2d[:, 3] ./= maximum(data_pa2d[:, 3])
const unc = 0.02

function find_parameters(data)
    ndata = size(data, 1)
    pas = [data[1, 1]]
    for i in 2:ndata
        pa = data[i, 1]
        if pa == pas[1]
            break
        end
        @assert(pa > pas[end])
        push!(pas, pa)
    end
    npas = length(pas)
    @assert(ndata % npas == 0)
    ndp = ndata ÷ npas
    dps = zeros(ndp)
    for i in 1:ndp
        dps[i] = data[(i - 1) * npas + 1, 2]
        for j in 1:npas
            k = (i - 1) * npas + j
            @assert(data[k, 2] == dps[i])
            @assert(data[k, 1] == pas[j])
        end
    end
    return pas, dps
end

const pas, dps = find_parameters(data_pa2d)
const npa, ndp = length(pas), length(dps)

const powers_2d = reshape(data_pa2d[:, 3], npa, ndp)
const sppowers_2d = powers_2d .- powers_2d[1:1, :] # Single pass power

const prefix = joinpath(@__DIR__, "imgs", "cs_tweezer_20200505")

function real_model(pa, dp, p, clip=true)
    # 0th order power:
    #   P0 = a3 * (1 - d3 * sin(b3 * atan(c3 * DP)^2))^2 * sin(b1 * sin(c1 * PA)^2)
    # 1st order power:
    #   P1 = a2 * sin(b2 * atan(c2 * DP)^2)^2
    b1, c1, a2, b2, c2, a3, b3, c3, d3 = p
    p0 = a3 * (1 - d3 * sin(b3 * atan(c3 * dp)^2))^2 * sin(b1 * sin(c1 * pa)^2)
    p1 = a2 * sin(b2 * atan(c2 * dp)^2)^2
    return clip ? min(p0 + p1, 1.0) : p0 + p1 # saturation
end

function model(x, p)
    function wrapper(x)
        pa = data_pa2d[x, 1]
        dp = data_pa2d[x, 2]
        return real_model(pa, dp, p)
    end
    return wrapper.(x)
end

function model2(x, p)
    b1, c1, a2, b2, c2, a3, b3, c3, d3 = p
    return a2 .* sin.(b2 .* atan.(c2 .* x).^2).^2
end
function model3(x, p, pa)
    b1, c1, a2, b2, c2, a3, b3, c3, d3 = p
    return a3 .* (1 .- d3 .* sin.(b3 .* atan.(c3 .* x).^2)).^2 * sin(b1 * sin(c1 * pa)^2)
end

function find_crossing(f1, f2, maxx)
    diff12(x) = f1(x) - f2(x)
    a0 = 0.0
    a1 = Float64(maxx)
    v0 = diff12(a0)
    v1 = diff12(a1)
    @assert(v0 * v1 < 0)
    while abs(a0 - a1) > 0.000001
        a2 = (a1 * v0 - a0 * v1) / (v0 - v1)
        v2 = diff12(a2)
        if v2 == 0
            return a2
        end
        a0, a1 = a1, a2
        v0, v1 = v1, v2
    end
    return a1, (f1(a1) + f2(a1)) / 2
end

function find_min(f, xmin, xmax)
    w1 = (sqrt(5) - 1) / 2
    w2 = 1 - w1

    a0 = Float64(xmin)
    a3 = Float64(xmax)
    a1 = a0 * w1 + a3 * w2
    a2 = a0 * w2 + a3 * w1
    v0 = f(a0)
    v1 = f(a1)
    v2 = f(a2)
    v3 = f(a3)
    @assert(min(v1, v2) < max(v0, v3))
    while abs(a0 - a1) > 0.00001
        if v1 > v2
            a0, a1, a2, a3 = a1, a2, a1 * w2 + a3 * w1, a3
            v0, v1, v2, v3 = v1, v2, f(a2), v3
        else
            a0, a1, a2, a3 = a0, a0 * w1 + a2 * w2, a1, a2
            v0, v1, v2, v3 = v0, f(a1), v1, v2
        end
    end
    return a1 > a2 ? (a2, v2) : (a1, v2)
end

find_crossing_23(fit, maxamp, pa) =
    find_crossing(x->model2(x, fit.param), x->model3(x, fit.param, pa), maxamp)

find_min_23(fit, maxamp, pa) =
    find_min(x->model2(x, fit.param) + model3(x, fit.param, pa), 0, maxamp)

function find_max_pa(fit, maxamp)
    x, y = find_min(x->-real_model(x, 0, fit.param, false), 0, maxamp)
    return x, -y
end

fit_pa2d = fit_data(model, 1:size(data_pa2d, 1), data_pa2d[:, 3],
                    [1.8, 2.1, 0.8, 1.6, 2.4, 1.1, 2.3, 1.8, 0.76], plotx=false)

figure()
pcolormesh(reshape(data_pa2d[:, 1], length(pas), length(dps)),
           reshape(data_pa2d[:, 2], length(pas), length(dps)), powers_2d; shading="flat")
xlabel("PAAOM/AMP")
ylabel("PADPaom/AMP")
colorbar()
title("Total power (normalized)")
NaCsPlot.maybe_save("$(prefix)_2d")

figure()
for i in 1:ndp
    errorbar(pas, sppowers_2d[:, i], fill(unc, npa))
end
xlim([0, pas[end] * 1.04])
ylim([0, 1])
xlabel("PAAOM/AMP")
title("Single Pass Power (normalized)")
grid()
NaCsPlot.maybe_save("$(prefix)_sp")

figure()
for i in 1:ndp
    errorbar(pas[2:end], sppowers_2d[2:end, i] ./ sppowers_2d[2:end, 4], fill(unc, npa - 1))
end
xlim([pas[2] * 0.95, pas[end] * 1.04])
xlabel("PAAOM/AMP")
title("Single Pass Power Ratio")
gca()[:set_yscale]("log", nonposy="clip")
grid()
NaCsPlot.maybe_save("$(prefix)_sp_ratio")

param_str(v) = @sprintf("%.4f", v)

plot_dps = linspace(0, dps[end] * 1.02, 1000)

figure()
for i in 1:npa
    errorbar(dps, powers_2d[i, :], fill(unc, ndp), fmt="C$i.", label="$(pas[i])")
    plot(plot_dps, real_model.(pas[i], plot_dps, (fit_pa2d.param,)), "C$i")
end
xlim([0, dps[end] * 1.04])
ylim([0, 1])
legend(fontsize="x-small", ncol=4, borderpad=0.2, labelspacing=0.2,
       handletextpad=0.3, columnspacing=0.2, borderaxespad=0.4, loc="upper right")
b1, c1, a2, b2, c2, a3, b3, c3, d3 = fit_pa2d.param
text(0.73, 0.50, "\$b_1=$(param_str(b1))\$", fontsize="small", color="C0")
text(0.73, 0.43, "\$c_1=$(param_str(c1))\$", fontsize="small", color="C0")
text(0.73, 0.30, "\$a_2=$(param_str(a2))\$", fontsize="small", color="C0")
text(0.40, 0.23, "\$b_2=$(param_str(b2))\$", fontsize="small", color="C0")
text(0.73, 0.23, "\$c_2=$(param_str(c2))\$", fontsize="small", color="C0")
text(0.40, 0.10, "\$a_3=$(param_str(a3))\$", fontsize="small", color="C0")
text(0.73, 0.10, "\$b_3=$(param_str(b3))\$", fontsize="small", color="C0")
text(0.40, 0.03, "\$c_3=$(param_str(c3))\$", fontsize="small", color="C0")
text(0.73, 0.03, "\$d_3=$(param_str(d3))\$", fontsize="small", color="C0")
xlabel("PADPaom/AMP")
title("Total Power (normalized)")
grid()
NaCsPlot.maybe_save("$(prefix)_total")

plot_pas = linspace(0, pas[end] * 1.02, 1000)

xmax_pa, ymax_pa = find_max_pa(fit_pa2d, 1)
figure()
errorbar(pas, powers_2d[:, 1], fill(unc, npa), fmt="C0.")
plot(plot_pas, real_model.(plot_pas, 0, (fit_pa2d.param,), false), "C0")
xlim([0, pas[end] * 1.04])
ylim([0, 1.4])
plot(xmax_pa, ymax_pa, "rs")
ax = gca()
ax[:annotate]("\$x_{max}=$(@sprintf("%.5f", xmax_pa))\$\n\$y_{max}=$(@sprintf("%.4f", ymax_pa))\$",
              (xmax_pa, ymax_pa + 0.004), xytext=(0.4, 0.6),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"r"),
              color="r", fontsize="small")
xlabel("PAAOM/AMP")
title("PAAOM")
grid()
NaCsPlot.maybe_save("$(prefix)_paaom")

x0, y0 = find_crossing_23(fit_pa2d, 0.6, xmax_pa)
xmin, ymin = find_min_23(fit_pa2d, 0.6, xmax_pa)
figure()
plot(plot_dps, model2(plot_dps, fit_pa2d.param), "C0-", label="1st order")
plot(plot_dps, model3(plot_dps, fit_pa2d.param, xmax_pa), "C1-", label="0th order")
plot(plot_dps, model2(plot_dps, fit_pa2d.param) .+ model3(plot_dps, fit_pa2d.param, xmax_pa),
     "--", color="#d452d4", label="Both")
plot(x0, y0, "bs")
plot(xmin, ymin, "rs")
ax = gca()
ax[:annotate]("x=$(@sprintf("%.5f", x0))\ny=$(@sprintf("%.4f", y0))",
              (x0 + 0.004, y0), xytext=(0.55, 0.12),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"b"),
              color="b", fontsize="small")
ax[:annotate]("\$x_{min}=$(@sprintf("%.5f", xmin))\$\n\$y_{min}=$(@sprintf("%.4f", ymin))\$",
              (xmin, ymin + 0.004), xytext=(0.235, 0.82),
              arrowprops=Dict(:arrowstyle=>"fancy", :color=>"r"),
              color="r", fontsize="small")
title("PADPaom")
legend(ncol=1, fontsize="x-small")
grid()
xlim([0, dps[end] * 1.05])
ylim([0, 1.3])
xlabel("PADPaom/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_padpaom")

NaCsPlot.maybe_show()
