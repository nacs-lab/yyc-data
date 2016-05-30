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
            ψ = exp(-((i * grid_space - x_center + 0.25f0) / 0.2f0)^2)
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

    function run(param)
        β, det_hz, trapf = param
        # β = 0.6
        totalt = 8_000f0
        println("Running β=$β, δ=$det_hz, trapf=$trapf")
        # m here is actually m / ħ
        # m_Na = Float32(22.98977e-3 / 6.02214129e23 /
        #                (1.0545717253362894e-34 * 1e6))
        m_Cs = Float32(132.905451933e-3 / 6.02214129e23 /
                       (1.0545717253362894e-34 * 1e6))
        ω_g = Float32(2π * trapf) # f = 160kHz
        ω_e = sign(β) * ω_g * √(abs(β))
        h_trap = HTrap{Float32}(m_Cs, (ω_g, ω_e))

        λ_res = Float32(0.852)

        # k, Γ
        o_decay = OpticalDecay{Float32}(2π / λ_res, 2π * 5)

        trap_depth = 1.2 * 2π * 20 * (trapf / 0.16)^2 # 2mK
        e_trap_depth = trap_depth * ω_e^2 / ω_g^2
        δ_offset = trap_depth - e_trap_depth

        δ_0 = 2π * det_hz

        # k, Ω, δ, τ_θ
        δ = δ_0 - δ_offset
        Ω = 2π * 5 / 5
        o_drive1 = OpticalDrive{Float32}(2π / λ_res, Ω, δ, 1000.0)
        o_drive2 = OpticalDrive{Float32}(-2π / λ_res, Ω, δ, 1000.0)
        # o_drive3 = OpticalDrive{Float32}(0, Ω * 2, δ, 1000.0)

        h_system = HSystem(h_trap, o_decay, (o_drive1, o_drive2))
        # h_system = HSystem(h_trap, o_decay, (o_drive1,))
        # h_system = HSystem(h_trap, o_decay, (o_drive2,))

        grid_size = 512
        grid_space = 0.0025f0
        tstep = 0.001f0
        nstep = round(Int, totalt ÷ tstep)
        p_sys = SystemPropagator(h_system, tstep, grid_space,
                                 nstep, grid_size)
        e_thresh = trap_depth
        ψ0 = gen_ψ0(grid_size, grid_space)
        _accum = EnergyRecorder(p_sys, e_thresh)
        monte_carlo = 50
        accum = MonteCarloAccumulator(_accum, monte_carlo)
        propagate(p_sys, ψ0, accum)
        return accum
    end
end

println("start")

# βs = [1.67, 1, 0.6, -0.6]
# βs = [2.5 ,1]
# βs = [4, 2.5]
βs = [1, 0.6]
det_list = [
            # [-90, -85, -80, -75, -70, -65, -60, -55, -50, -45],
            # [-40, -38.75, -37.5, -36.25, -35, -30],
            # [-20, -18.75, -17.5, -17.1875, -16.875, -16.5625, -16.25, -15],
            [-10, -5, -2.5, -2.1875, -1.875, -1.5625, -1.25, 0, 5, 10],
            [-20, -15, -10, -5, -2.5, 0, 2.5, 5, 10, 15],
            # [-35, -30, -25, -20, -15, -10, -5, 0,
            #  5, 10, 15, 20, 25, 30, 35, 40],
            ]
trapfs = [0.32, 0.16, 0.08]
# ts = [10000, 6000, 3000]
# ts = [8000]

# t = 10_000, β = 1
# params = [-15, -10, -7.5, -6.25, -5, 0, 5, 10]
# t = 10_000, β = 0.6
# params = [-20, -15, -10, -5, 0, 5, 10, 15]
# t = 6_000, β = 1
# params = [-15, -10, -7.5, -6.25, -5, 0, 5, 10]
# t = 6_000, β = 0.6
# params = [-20, -15, -10, -5, 0, 5, 10, 15]
# t = 6_000, β = 1.67
# params = [-30, -25, -20, -17.5, -15, -10, -5, 0]
# t = 6_000, β = -0.6
# params = [-35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40]
params = Vector{NTuple{3,Float32}}()
for trapf in trapfs
    for i in 1:length(βs)
        β = βs[i]
        for det in det_list[i]
            push!(params, (β, det, trapf))
        end
    end
end

xax_name = "Free space detuning"
@time accums = pmap(run, params)

# exit()

using PyPlot

# for _accum in accums
#     # println(_accum)
#     plot_accum(_accum)
# end
# ylim(0, ylim()[2] * 1.1)

final_e = Float32[_accum.Es[end] for _accum in accums]
final_e2 = Float32[_accum.Es2[end] for _accum in accums]
println(("final_energy", params, final_e, final_e2))
# figure()
# errorbar(params, final_e, final_e2)
# xlabel(xax_name)
# ylabel("Final energy")
# grid()
# ylim(0, ylim()[2])
# savefig("final_energy.png")
# close()

t_esc = Float32[_accum.t_esc for _accum in accums]
t_esc2 = Float32[_accum.t_esc2 for _accum in accums]
println(("escape_time", params, t_esc, t_esc2))
# figure()
# errorbar(params, t_esc, t_esc2)
# xlabel(xax_name)
# ylabel("Escape time (us)")
# grid()
# ylim(0, ylim()[2])
# savefig("escape_time.png")
# close()

pcount = Float32[_accum.pcount for _accum in accums]
pcount2 = Float32[_accum.pcount2 for _accum in accums]
println(("photon_count", params, pcount, pcount2))
# figure()
# errorbar(params, pcount, pcount2)
# xlabel(xax_name)
# ylabel("Photon count")
# grid()
# ylim(0, ylim()[2])
# savefig("photon_count.png")
# close()

# show()
