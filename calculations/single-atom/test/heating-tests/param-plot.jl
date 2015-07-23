#!/usr/bin/julia -f

@everywhere begin
    # include(joinpath(dirname(@__FILE__), "two-level-1d.jl"))

    # Time unit: μs
    # Length unit: μm
    # Frequency unit: MHz
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

    function run(ratio)
        println("Running $ratio")
        # m here is actually m / ħ
        m_Na = Float32(22.98977e-3 / 6.02214129e23 /
                       (1.0545717253362894e-34 * 1e6))
        ω_g = Float32(2π * 0.2) # f = 200kHz
        ω_e = Float32(sign(ratio) * ω_g * sqrt(abs(ratio)))
        h_trap = HTrap{Float32}(m_Na, (ω_g, ω_e))

        λ_res = Float32(0.589)

        # k, Γ
        o_decay = OpticalDecay{Float32}(2π / λ_res, 2π * 10.0)

        # k, Ω, δ, τ_θ
        δ = -2π * 5.0 * 6
        Ω = 2π * 2.5
        o_drive1 = OpticalDrive{Float32}(2π / λ_res, Ω, δ, 1000.0)
        o_drive2 = OpticalDrive{Float32}(-2π / λ_res, Ω, δ, 1000.0)

        h_system = HSystem(h_trap, o_decay, (o_drive1, o_drive2))

        grid_size = 512
        grid_space = 0.005f0
        p_sys = SystemPropagator(h_system, 0.0025f0, grid_space,
                                 80000, grid_size)
        ψ0 = gen_ψ0(grid_size, grid_space)
        _accum = EnergyRecorder(p_sys)
        monte_carlo = 100
        accum = MonteCarloAccumulator(_accum, monte_carlo)
        propagate(p_sys, ψ0, accum)
        return accum
    end
end

println("start")

ratios = linspace(-5f0, 10f0, 8)
@time accums = pmap(run, ratios)

using PyPlot

for _accum in accums
    plot_accum(_accum)
end

figure()
final_e = Float32[_accum.Es[end] for _accum in accums]
final_e2 = Float32[_accum.Es2[end] for _accum in accums]
errorbar(ratios, final_e, final_e2)

show()
