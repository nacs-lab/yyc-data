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

    function run(ratio)
        println("Running $ratio")
        # m here is actually m / ħ
        m_Na = Float32(22.98977e-3 / 6.02214129e23 /
                       (1.0545717253362894e-34 * 1e6))
        ω_g = Float32(2π * 0.2) # f = 200kHz
        ω_e = ω_g
        h_trap = HTrap{Float32}(m_Na, (ω_g, ω_e))

        λ_res = Float32(0.589)

        # k, Γ
        o_decay = OpticalDecay{Float32}(2π / λ_res, 2π * 10.0)

        # k, Ω, δ, τ_θ
        δ = -2π * 0.0
        Ω = 2π * ratio
        o_drive1 = OpticalDrive{Float32}(2π / λ_res, Ω, δ, 1000.0)
        o_drive2 = OpticalDrive{Float32}(-2π / λ_res, Ω, δ, 1000.0)

        h_system = HSystem(h_trap, o_decay, (o_drive1, o_drive2))
        h_system = HSystem(h_trap, o_decay, (o_drive1,))
        # h_system = HSystem(h_trap, o_decay, (o_drive2,))

        grid_size = 512
        grid_space = 0.005f0
        p_sys = SystemPropagator(h_system, 0.005f0, grid_space,
                                 600000, grid_size)
        e_thresh = (maximum(p_sys.E_k) + maximum(p_sys.E_x[1])) / 4
        ψ0 = gen_ψ0(grid_size, grid_space)
        _accum = EnergyRecorder(p_sys, e_thresh)
        monte_carlo = 100
        accum = MonteCarloAccumulator(_accum, monte_carlo)
        propagate(p_sys, ψ0, accum)
        return accum
    end
end

println("start")

ratios = linspace(0.5, 2, 4) * 1.25
@time accums = pmap(run, ratios)

using PyPlot

for _accum in accums
    plot_accum(_accum)
end

figure()
final_e = Float32[_accum.Es[end] for _accum in accums]
final_e2 = Float32[_accum.Es2[end] for _accum in accums]
errorbar(ratios, final_e, final_e2)
xlabel("Rabi frequency")
ylabel("Final energy")
grid()
ylim(0, ylim()[2])

figure()
t_esc = Float32[_accum.t_esc for _accum in accums]
t_esc2 = Float32[_accum.t_esc2 for _accum in accums]
errorbar(ratios, t_esc, t_esc2)
xlabel("Rabi frequency")
ylabel("Escape time")
grid()
ylim(0, ylim()[2])

figure()
pcount = Float32[_accum.pcount for _accum in accums]
pcount2 = Float32[_accum.pcount2 for _accum in accums]
errorbar(ratios, pcount, pcount2)
xlabel("Rabi frequency")
ylabel("Photon count")
grid()
ylim(0, ylim()[2])

show()
