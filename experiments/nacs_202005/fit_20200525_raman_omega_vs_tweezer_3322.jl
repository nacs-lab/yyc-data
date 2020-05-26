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

const freqs = [560 15 4.023 0.059
               560 6 1.182 0.033
               560 3 0.488 0.015

               503 15 2.073 0.051
               503 9 1.074 0.036
               503 6 0.680 0.022
               503 3 0.2544 0.0077

               605 15 6.39 0.15
               605 6 1.985 0.043
               605 3 0.845 0.042
               ]

function model_omega(x, p)
    return p[1] .* x.^1.29
end

broacast_array(f, x) = f(x)
broacast_array(f, x::AbstractArray) = f.(x)

# x: freq, power
# p: offset, strength
function gen_model(fpa0)
    function model(xs, p)
        function real_model(x)
            freq, power = x
            offset, strength = p
            p_omega = (offset - strength / (freq - fpa0),)
            return model_omega(power, p_omega)
        end
        return broacast_array(real_model, xs)
    end
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

# function fit_freq_model(model, data, _p0)
#     p0 = [_p0; zeros(length(data.freqs))]
#     return fit_data(model, data.x, data.y, data.unc, p0, plotx=false)
# end

# function get_plot_data_freq(data, fit, model, freq)
#     local freq_id
#     xs = Float64[]
#     ys = Float64[]
#     uncs = Float64[]
#     for i in 1:length(data.x)
#         x = data.x[i]
#         x[1] == freq || continue
#         if !@isdefined(freq_id)
#             freq_id = x[3]
#         else
#             @assert(freq_id == x[3])
#         end
#         push!(xs, x[2])
#         push!(ys, data.y[i])
#         push!(uncs, data.unc[i])
#     end
#     function plot_func(x)
#         return model((freq, x, freq_id), fit.param)
#     end
#     plotx = get_plot_range(xs)
#     return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
# end

function get_plot_data_power(data, fit, model, power)
    xs = Float64[]
    xids = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 1:length(data.x)
        x = data.x[i]
        x[2] == power || continue
        push!(xs, x[1])
        push!(ys, data.y[i])
        push!(uncs, data.unc[i])
    end
    function plot_func(x)
        return model((x, power,), fit.param)
    end
    plotx = get_plot_range(data.freqs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

const model = gen_model(705)
const data = gen_data(freqs)
const fit1 = fit_data(model, data.x, data.y, data.unc, [0.1, 0.1], plotx=false)
@show fit1.uncs

const prefix = joinpath(@__DIR__, "imgs", "fit_20200525_raman_omega_vs_tweezer_3322")

figure()
const powers_plot = [15, 9, 6, 3]
# const plot_freqs = get_plot_range(data.freqs)
for i in 1:length(powers_plot)
    power = powers_plot[i]
    pd = get_plot_data_power(data, fit1, model, power)
    plot(pd.plotx, pd.ploty, "C$(i - 1)", label="$(power) mW")
    errorbar(pd.x, pd.y, pd.unc, fmt="C$(i - 1).")
end
# p: offset, strength
text(500, 5.0, ("\$\\left(a-\\dfrac{b}{f-705 GHz}\\right)\\cdot P^{1.29}\$"))
text(500, 2.8, ("\$a=$(fit1.uncs[1] * 1000)\$ Hz/mW\$^{1.29}\$\n" *
              "\$b=$(fit1.uncs[2])\$ kHz\$\\cdot\$GHz/mW\$^{1.29}\$"), fontsize="small")
# axhline((fit1.param[1] .- 770) .* 1000, color="c3", ls="--")
legend(ncol=2, fontsize="x-small")
grid()
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("\$\\Omega_{Raman} (2\\pi\\cdot \\mathrm{kHz})\$")
NaCsPlot.maybe_save("$(prefix)")

figure()
plotxf1 = linspace(400, 680, 1501)
plotxf2 = linspace(730, 1000, 1501)
plot(plotxf1, abs.(model.(tuple.(plotxf1, 15,), Ref(fit1.param))), "C0")
plot(plotxf2, abs.(model.(tuple.(plotxf2, 15,), Ref(fit1.param))), "C0")
grid()
title("15 mW")
ylim([0, 17])
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("\$\\Omega_{Raman} (2\\pi\\cdot \\mathrm{kHz})\$")
NaCsPlot.maybe_save("$(prefix)_2sides")

NaCsPlot.maybe_show()
