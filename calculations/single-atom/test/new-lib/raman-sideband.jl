#!/usr/bin/julia -f

using SingleAtom

# Time unit: μs
# Length unit: μm
# Frequency unit: MHz

# m here is actually m / ħ
m_Na = Float32(22.98977e-3 / 6.02214129e23 / (1.0545717253362894e-34 * 1e6))
ω_g = Float32(2π * 0.2) # f = 200kHz
ω_e = Float32(2π * 0.2) # f = 200kHz
λ_res = 0.589f0
k_res = Float32(2π / λ_res)
Γ_Na = Float32(2π * 10f0)

builder = SystemBuilder{Float32}()
add_state!(builder, :G1, :G1, 0)
add_state!(builder, :G2, :G2, 0)
add_state!(builder, :E, :E, 0)
add_potential!(builder, HarmonicPotential(ω_g), :G1)
add_potential!(builder, HarmonicPotential(ω_g), :G2)
add_potential!(builder, HarmonicPotential(ω_e), :E)

add_transition!(builder, :G1, :E, Transition{Trans_σ⁺}(k_res, Γ_Na, 1f0))
add_transition!(builder, :G2, :E, Transition{Trans_σ⁺}(k_res, Γ_Na / 2, 1f0))
add_transition!(builder, :G1, :G2, Transition{Trans_σ⁺}(k_res, 0f0, 1f0))

Ω = 2π * 0.5

add_drive!(builder, Drive{Vec3D{Complex64}(0, 2π * 0.05, 0)}(k_res, -ω_g, NaN,
                                                              1000.0),
           :G1, :G2)

add_drive!(builder, Drive{Vec3D{Complex64}(0, Ω, 0)}(-k_res, 0, NaN, 1000.0),
           :G2, :E)

# add_drive!(builder, Drive{Vec3D{Complex64}(0, 0.3Ω, 0)}(-k_res, 0, NaN, 1000.0),
#            :G1, :E)

# add_drive!(builder, Drive{Vec3D{Complex64}(0, Ω, 0)}(0, 0, NaN, 1000.0),
#            :G2, :E)

# add_drive!(builder, Drive{Vec3D{Complex64}(0, 0.3Ω, 0)}(0, 0, NaN, 1000.0),
#            :G1, :E)

sys = MotionSystem(Vec3D(1f0, 0f0, 0f0), m_Na, builder)

grid_size = 512
grid_space = 0.005f0

P = SystemPropagator(sys, 0.02f0, grid_space, 100000, grid_size)

e_thresh = (maximum(P.motion.E_k) + maximum(P.motion.E_x[1])) / 4

function gen_ψ0(grid_size, grid_space)
    x_center = (grid_size + 1) * grid_space / 2
    # Three states only
    ψ0 = StructOfArrays(Complex64, grid_size, 3)
    sum = 0.0
    @inbounds for i in 1:grid_size
        ψ = exp(-((i * grid_space - x_center + 0.3f0) / 0.2f0)^2)
        sum += abs2(ψ)
        ψ0[i, 1] = ψ
        ψ0[i, 2] = 0
        ψ0[i, 3] = 0
    end
    sum = sqrt(sum)
    @inbounds @simd for i in 1:grid_size
        ψ0[i, 1] /= sum
    end
    ψ0
end

ψ0 = gen_ψ0(grid_size, grid_space)

@enum PlotType PlotWFX PlotWFK PlotE

# const plot_type = PlotWFX
# const plot_type = PlotWFK
const plot_type = PlotE
const monte_carlo = 100

_measure = if plot_type == PlotWFX
    WaveFuncMeasure{SnapshotX}(P)
elseif plot_type == PlotWFK
    WaveFuncMeasure{SnapshotK}(P)
elseif plot_type == PlotE
    EnergyMeasure(P, :G1, e_thresh)
end

if monte_carlo > 1
    _measure = MonteCarloMeasure(_measure, monte_carlo)
end

setup = StaticSetup(ψ0)

println("start")

@time precompile(propagate, tuple(Base.typesof(P, setup,
                                               _measure).parameters...))
gc()
# @profile propagate(P, setup, _measure)
# Profile.print(C=true)
@time propagate(P, setup, _measure)

# exit()

plot_measure(_measure)

show()
