#!/usr/bin/julia -f

include("two-level-1d.jl")

# Time unit: μs
# Length unit: μm
# Frequency unit: MHz

# m here is actually m / ħ
m_Na = 22.98977e-3 / 6.02214129e23 / (1.0545717253362894e-34 * 1e6)
ω_g = 2π * 0.1 # f = 100kHz
ω_e = 2π * 0.1 # f = 100kHz
h_trap = HTrap(m_Na, (ω_g, ω_e))

# k, Γ
o_decay = OpticalDecay(2π * 0.589, 2π * 10.0)

# k, Ω, δ, τ_θ
o_drive1 = OpticalDrive(2π * 0.589, 2π * 5.0, 2π * 0.0, 100.0)

h_system = HSystem(h_trap, o_decay, (o_drive1,))

grid_size = 256
grid_space = 0.01
p_sys = SystemPropagator(h_system, 0.005, grid_space, 10000, grid_size)


function gen_ψ0(grid_size, grid_space)
    x_center = (grid_size + 1) * grid_space / 2
    ψ0 = Array{Complex128}(2, grid_size)
    sum = 0.0
    @inbounds for i in 1:grid_size
        ψ = exp(-((i * grid_space - x_center + 0.2) / 0.2)^2)
        sum += abs2(ψ)
        ψ0[1, i] = ψ
        ψ0[2, i] = 0
    end
    sum = sqrt(sum)
    @inbounds for i in 1:grid_size
        ψ0[1, i] /= sum
    end
    ψ0
end

ψ0 = gen_ψ0(grid_size, grid_space)

@enum PlotType PlotWFX PlotWFK PlotE

# const plot_type = PlotWFX
# const plot_type = PlotWFK
const plot_type = PlotE
const monte_carlo = 100

_accum = if plot_type == PlotWFX
    WaveFuncRecorder{AccumX}(p_sys)
elseif plot_type == PlotWFK
    WaveFuncRecorder{AccumK}(p_sys)
elseif plot_type == PlotE
    EnergyRecorder(p_sys)
end

if monte_carlo > 1
    _accum = MonteCarloAccumulator(_accum, monte_carlo)
end

println("start")

precompile(propagate, tuple(Base.typesof(p_sys, ψ0, _accum).parameters...))
# @time propagate(p_sys, ψ0, _accum)
gc()
@time propagate(p_sys, ψ0, _accum)

using PyPlot

function plot_accum_img(img::Matrix{Float64})
    xsize, ysize = size(img)

    if xsize > ysize * 3
        xscale = xsize ÷ (ysize * 2)
        img = img[1:xscale:end, :]
    elseif ysize > xsize * 3
        yscale = ysize ÷ (xsize * 2)
        img = img[:, 1:yscale:end]
    end

    figure()
    imshow(img)
    colorbar()

    figure()
    imshow(log(img))
    colorbar()

    nothing
end

function plot_accum(accum::WaveFuncRecorder)
    ψs = accum.ψs

    img = Array{Float64}(size(ψs, 2, 3))

    for i in 1:size(img, 2)
        @inbounds for j in 1:size(img, 1)
            img[j, i] = abs2(ψs[1, j, i]) + abs2(ψs[2, j, i])
        end
    end
    plot_accum_img(img)
end

function plot_accum(accum::WaveFuncMonteCarloRecorder)
    ψs = accum.ψs2

    img = Array{Float64}(size(ψs, 2, 3))

    for i in 1:size(img, 2)
        @inbounds for j in 1:size(img, 1)
            img[j, i] = ψs[1, j, i] + ψs[2, j, i]
        end
    end
    plot_accum_img(img)
end

function plot_accum(accum::EnergyRecorder)
    figure()
    plot(accum.Es)
    grid()
    ylim(0, ylim()[2] * 1.1)
end

function plot_accum(accum::EnergyMonteCarloRecorder)
    figure()
    errorbar(1:length(accum.Es), accum.Es, accum.Es2)
    grid()
    ylim(0, ylim()[2] * 1.1)
end

plot_accum(_accum)
show()
