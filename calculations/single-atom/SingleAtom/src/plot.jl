#!/usr/bin/julia -f

module Plotting

using PyPlot
using ..Measure

export plot_measure

function plot_measure_img(img::Matrix{Float64})
    xsize, ysize = size(img)

    if xsize > ysize * 3
        xscale = xsize ÷ (ysize * 2)
        img = img[1:xscale:end, :]
    elseif ysize > xsize * 3
        yscale = ysize ÷ (xsize * 2)
        img = img[:, 1:yscale:end]
    end

    imshow(img)
    colorbar()

    nothing
end

function plot_measure(measure::WaveFuncMeasure)
    ψs = measure.ψs

    nstates, = size(ψs, 2)
    img = zeros(Float64, size(ψs, 1, 3))

    @inbounds for i in 1:size(img, 2)
        for k in 1:nstates
            @simd for j in 1:size(img, 1)
                img[j, i] += abs2(ψs[j, k, i])
            end
        end
    end
    plot_measure_img(img)
end

function plot_measure(measure::WaveFuncMonteCarloMeasure)
    ψs = measure.ψs2

    nstates, = size(ψs, 2)
    img = zeros(Float64, size(ψs, 1, 3))

    @inbounds for i in 1:size(img, 2)
        for k in 1:nstates
            @simd for j in 1:size(img, 1)
                img[j, i] += ψs[j, k, i]
            end
        end
    end
    plot_measure_img(img)
end

function plot_measure(measure::EnergyMeasure)
    plot(measure.Es)
    grid()
    ylim(0, ylim()[2] * 1.1)
end

function plot_measure(measure::EnergyMonteCarloMeasure)
    @printf("Escape time: %.2f±%.2f\n", measure.t_esc.v, measure.t_esc.v2)
    @printf("Photon Emitted: %.2f±%.2f\n", measure.pcount.v,
            measure.pcount.v2)
    errorbar((1:length(measure.Es)) * measure.sub_measure.dt,
             measure.Es, measure.Es2)
    axvline(x=measure.t_esc.v, color="r")
    axvline(x=measure.t_esc.v - measure.t_esc.v2, color="b", linestyle="--")
    axvline(x=measure.t_esc.v + measure.t_esc.v2, color="b", linestyle="--")
    grid()
    ylim(0, ylim()[2] * 1.1)
end

end
