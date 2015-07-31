#!/usr/bin/julia -f

include("two-level-1d.jl")

# Time unit: μs
# Length unit: μm
# Frequency unit: MHz

# m here is actually m / ħ
m_Na = Float32(22.98977e-3 / 6.02214129e23 / (1.0545717253362894e-34 * 1e6))
ω_g = Float32(2π * 0.2) # f = 200kHz
ω_e = Float32(2π * 0.2) # f = 200kHz
h_trap = HTrap(m_Na, (ω_g, ω_e))

λ_res = Float32(0.589)

# k, Γ
o_decay = OpticalDecay{Float32}(2π / λ_res, 2π * 10.0)

# k, Ω, δ, τ_θ
δ = -2π * 0.0
Ω = 2π * 1.25
o_drive1 = OpticalDrive{Float32}(2π / λ_res, Ω, δ, 1000.0)
o_drive2 = OpticalDrive{Float32}(-2π / λ_res, Ω, δ, 1000.0)

# h_system = HSystem(h_trap, o_decay, (o_drive1, o_drive2))
h_system = HSystem(h_trap, o_decay, (o_drive1,))
# h_system = HSystem(h_trap, o_decay, (o_drive2,))

grid_size = 512
grid_space = 0.005f0
p_sys = SystemPropagator(h_system, 0.005f0, grid_space, 100000, grid_size)
e_thresh = (maximum(p_sys.E_k) + maximum(p_sys.E_x[1])) / 4

plot_title = "\$\\omega_e = $ω_e\$, \$\\omega_g = $ω_g\$"

function gen_ψ0(grid_size, grid_space)
    x_center = (grid_size + 1) * grid_space / 2
    ψ0 = Array{Complex64}(2, grid_size)
    sum = 0.0
    @inbounds for i in 1:grid_size
        ψ = exp(-((i * grid_space - x_center + 0.5f0 * 0) / 0.2f0)^2)
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
    EnergyRecorder(p_sys, e_thresh)
end

if monte_carlo > 1
    _accum = MonteCarloAccumulator(_accum, monte_carlo)
end

println("start")

@time precompile(propagate, tuple(Base.typesof(p_sys, ψ0, _accum).parameters...))
# @time propagate(p_sys, ψ0, _accum)
# @code_native propagate(p_sys, ψ0, _accum)
gc()
@time propagate(p_sys, ψ0, _accum)

using PyPlot

plot_accum(_accum)
title(plot_title)

show()
