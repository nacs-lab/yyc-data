#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival, get_plot_range
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const freqs = [668.35 21 772.1367 0.0082
               668.35 21 772.1376 0.0049
               668.35 21 772.1172 0.0048
               668.35 18 771.924 0.033
               668.35 12 771.323 0.050

               560 15 770.59417 0.00013
               560 6 770.369081 0.000037
               560 3 770.287276 0.000015

               565 15 770.61656 0.00021
               565 6 770.378274 0.000080
               565 3 770.291891 0.000038

               555 15 770.57744 0.00013
               555 6 770.362165 0.000043
               555 3 770.283665 0.000017

               503 15 770.483553 0.000049
               503 12 770.43248 0.00011
               503 9 770.379505 0.000028
               503 6 770.324463 0.000019
               503 3 770.2648307 0.0000084

               492 15 770.47571 0.00010
               492 12 770.42599 0.00032
               492 9 770.375007 0.000069
               492 6 770.321355 0.000037
               492 3 770.263402 0.000028

               510 15 770.50117 0.00083
               510 12 770.44659 0.00021
               510 9 770.390205 0.000079
               510 6 770.331510 0.000057
               510 3 770.258511 0.000023

               605 15 770.76208 0.00027
               605 12 770.65624 0.00046
               605 9 770.54816 0.00058
               605 6 770.437052 0.000072
               605 3 770.321367 0.000038
               ]

function model_sqr(x, p)
    return p[1] .+ (p[2] .+ p[3] .* x) .* x
end

broacast_array(f, x) = f(x)
broacast_array(f, x::AbstractArray) = f.(x)

# x: freq, power, idx
# p: framan0, fpa0, offset, strength, bs...
function model1(xs, p)
    function real_model(x)
        freq, power, idx = x
        framan0, fpa0, offset, strength = p
        b = idx == 0 ? 0.0 : p[4 + idx]
        p_sqr = (framan0, offset - strength / (freq - fpa0), b)
        return model_sqr(power, p_sqr)
    end
    return broacast_array(real_model, xs)
end

function gen_data(data_in)
    freq_ids = Dict{Float64,Int}()
    xs = Tuple{Float64,Float64,Int}[]
    ys = Float64[]
    uncs = Float64[]
    idmax = 0
    for i in 1:size(data_in, 1)
        freq = data_in[i, 1]
        power = data_in[i, 2]
        freq_id = get!(freq_ids, freq) do
            idmax += 1
            return idmax
        end
        push!(xs, (freq, power, freq_id))
        push!(ys, data_in[i, 3])
        push!(uncs, data_in[i, 4])
    end
    return (x=xs, y=ys, unc=uncs, freqs=collect(keys(freq_ids)))
end

function fit_freq_model(model, data, _p0)
    p0 = [_p0; zeros(length(data.freqs))]
    return fit_data(model, data.x, data.y, data.unc, p0, plotx=false)
end

function get_plot_data_freq(data, fit, model, freq)
    local freq_id
    xs = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 1:length(data.x)
        x = data.x[i]
        x[1] == freq || continue
        if !@isdefined(freq_id)
            freq_id = x[3]
        else
            @assert(freq_id == x[3])
        end
        push!(xs, x[2])
        push!(ys, data.y[i])
        push!(uncs, data.unc[i])
    end
    function plot_func(x)
        return model((freq, x, freq_id), fit.param)
    end
    plotx = get_plot_range(xs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

# function get_plot_data_power(data, fit, model, power)
#     xs = Float64[]
#     xids = Float64[]
#     ys = Float64[]
#     uncs = Float64[]
#     for i in 1:length(data.x)
#         x = data.x[i]
#         x[2] == power || continue
#         push!(xs, x[1])
#         push!(ys, data.y[i])
#         push!(uncs, data.unc[i])
#     end
#     function plot_func(x)
#         return model((x, power, 0), fit.param)
#     end
#     plotx = get_plot_range(xs)
#     return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
# end

const data = gen_data(freqs)
const fit1 = fit_freq_model(model1, data, [770.202, 705, 0, 100])
@show fit1.uncs

const prefix = joinpath(@__DIR__, "imgs", "fit_20200508_raman_vs_tweezer_3322")

figure()
errorbar(data.freqs, fit1.param[5:end] * 1000, fit1.unc[5:end] * 1000, fmt="C0.")
grid()
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("Quadratic shift (\$\\mathrm{kHz/mW^2}\$)")
title("Quadratic shift")
NaCsPlot.maybe_save("$(prefix)_2")

const freqs_sorted = sort(data.freqs)
const nfreq = length(data.freqs)
const nplotx = ceil(Int, sqrt(nfreq))
const nploty = ceil(Int, nfreq / nplotx)

figure(figsize=[6.3 * nplotx, 5.6 * nploty])
for i in 1:nfreq
    freq = freqs_sorted[i]
    subplot(nploty, nplotx, i)
    plot1 = get_plot_data_freq(data, fit1, model1, freq)
    errorbar(plot1.x, (plot1.y .- 770) .* 1000, plot1.unc .* 1000, fmt="C0s")
    plot(plot1.plotx, (plot1.ploty .- 770) .* 1000, "C0")
    grid()
    xlabel("Tweezer Power (mW)")
    ylabel("Raman Resonance (770XXX kHz)")
    title("$(freq + 288000) GHz")
end

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_ps")

function plot_freq_func(model, fit, power, x)
    return model((x, power, 0), fit.param)
end
figure()
const powers_plot = [15, 6, 3]
const plot_freqs = get_plot_range(data.freqs)
for i in 1:length(powers_plot)
    power = powers_plot[i]
    plot(plot_freqs, (plot_freq_func.(model1, Ref(fit1), power, plot_freqs) .- 770) .* 1000,
         "C$(i - 1)", label="$(power) mW")
end
# p: framan0, fpa0, offset, strength, bs...
text(560, 1700, ("\$f_{Raman0} + a\\cdot P\$\n" *
                 "  \$-\\dfrac{b}{f-f_{PA0}}\\cdot P\$"))
text(480, 740, ("\$f_{Raman0}=$(fit1.uncs[1])\$ MHz\n" *
                "\$f_{PA0}=$(fit1.uncs[2] + 288000)\$ GHz\n" *
                "\$a=$(fit1.uncs[3] * 1000)\$ kHz/mW\n" *
                "\$b=$(fit1.uncs[4])\$ MHz\$\\cdot\$GHz/mW"), fontsize="small")
axhline((fit1.param[1] .- 770) .* 1000, color="C3", ls="--")
legend(fontsize="small")
grid()
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("Raman Resonance (770XXX kHz)\${}_{linear\\ shift}\$")
NaCsPlot.maybe_save("$(prefix)_fs1")

# figure()
# errorbar(tweezer, freq, freq_unc, fmt="C0s")
# plot(fit_freq.plotx, fit_freq.ploty, "C0", label="With 668.35 GHz point")
# plot(fit_freq2.plotx, fit_freq2.ploty, "C1", label="\$f_{res}=705\\ GHz\$")
# # plot(fit_freq3.plotx, fit_freq3.ploty, "C2", label="\$f_{res}=705, -15\\ GHz\$")
# legend(fontsize="small")
# grid()
# xlabel("Tweezer Frequency (288XXX GHz)")
# ylabel("Raman Resonance (770XXX kHz)")
# title("Raman Resonance (3 mW)")

NaCsPlot.maybe_show()
