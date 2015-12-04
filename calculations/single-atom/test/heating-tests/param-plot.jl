#!/usr/bin/julia -f

@everywhere begin
    # Time unit: μs
    # Length unit: μm
    # Frequency unit: MHz
    function gen_ψ0(grid_size, grid_space)
        x_center = (grid_size + 1) * grid_space / 2
        ψ0 = Array{Complex64}(2, grid_size)
        sum = 0.0
        @inbounds for i in 1:grid_size
            ψ = exp(-((i * grid_space - x_center + 0.8f0 * 0) / 0.2f0)^2)
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

    function run(det_hz)
        println("Running $det_hz")
        # m here is actually m / ħ
        # m_Na = Float32(22.98977e-3 / 6.02214129e23 /
        #                (1.0545717253362894e-34 * 1e6))
        m_Cs = Float32(132.905451933e-3 / 6.02214129e23 /
                       (1.0545717253362894e-34 * 1e6))
        ω_g = Float32(2π * 0.07) # f = 70kHz
        ω_e = ω_g * √(0.6)
        h_trap = HTrap{Float32}(m_Cs, (ω_g, ω_e))

        λ_res = Float32(0.852)

        # k, Γ
        o_decay = OpticalDecay{Float32}(2π / λ_res, 2π * 5.234)

        trap_depth = 0.7 * 2π * 20 # 0.7mK
        e_trap_depth = trap_depth * ω_e^2 / ω_g^2
        δ_offset = trap_depth - e_trap_depth

        δ_0 = 2π * det_hz

        # k, Ω, δ, τ_θ
        δ = δ_0 + δ_offset
        Ω = 2π * 5.234 * 0.05
        o_drive1 = OpticalDrive{Float32}(2π / λ_res, Ω, δ, 1000.0)
        o_drive2 = OpticalDrive{Float32}(-2π / λ_res, Ω, δ, 1000.0)

        h_system = HSystem(h_trap, o_decay, (o_drive1, o_drive2))
        # h_system = HSystem(h_trap, o_decay, (o_drive1,))
        # h_system = HSystem(h_trap, o_decay, (o_drive2,))

        grid_size = 512
        grid_space = 0.005f0
        p_sys = SystemPropagator(h_system, 0.01f0, grid_space,
                                 100_000, grid_size)
        e_thresh = trap_depth
        ψ0 = gen_ψ0(grid_size, grid_space)
        _accum = EnergyRecorder(p_sys, e_thresh)
        monte_carlo = 100
        accum = MonteCarloAccumulator(_accum, monte_carlo)
        propagate(p_sys, ψ0, accum)
        return accum
    end
end

println("start")

params = linspace(-25, 10, 8)
xax_name = "Free space detuning"
@time accums = pmap(run, params)

# exit()

using PyPlot

for _accum in accums
    plot_accum(_accum)
end

figure()
final_e = Float32[_accum.Es[end] for _accum in accums]
final_e2 = Float32[_accum.Es2[end] for _accum in accums]
errorbar(params, final_e, final_e2)
xlabel(xax_name)
ylabel("Final energy")
grid()
ylim(0, ylim()[2])

figure()
t_esc = Float32[_accum.t_esc for _accum in accums]
t_esc2 = Float32[_accum.t_esc2 for _accum in accums]
errorbar(params, t_esc, t_esc2)
xlabel(xax_name)
ylabel("Escape time (us)")
grid()
ylim(0, ylim()[2])

figure()
pcount = Float32[_accum.pcount for _accum in accums]
pcount2 = Float32[_accum.pcount2 for _accum in accums]
errorbar(params, pcount, pcount2)
xlabel(xax_name)
ylabel("Photon count")
grid()
ylim(0, ylim()[2])

show()
