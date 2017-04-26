#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 15,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

import NaCsSim: Setup, System
import NaCsCalc: Trap
import NaCsCalc.Utils: interactive

const BuilderT = Setup.SeqBuilder{System.StateC,Void}
const sz = 500, 100, 100
statec() = System.StateC(sz...)

const m_Na = 23f-3 / 6.02f23
const k_Na = Float32(2π) / 589f-9
η(freq) = Trap.η(m_Na, freq, k_Na)

const η_op = (η(67f3), η(420f3), η(580f3))
const η_op_dri = (0f0, η(420f3) / sqrt(2f0), η(580f3) / sqrt(2f0))

const η_raman1 = (η(67f3) / sqrt(2f0), η(420f3) / sqrt(2f0), η(580f3) / sqrt(2f0))
const η_raman2 = (0f0, η(420f3) * sqrt(2f0), 0f0)
const η_raman3 = (0f0, 0f0, η(580f3) * sqrt(2f0))
const η_ramans = (η_raman1, η_raman2, η_raman3)

function threadmap(f, arg)
    # A really simple work queue
    n = length(arg)
    counter = Threads.Atomic{Int}(1)
    T = Core.Inference.return_type(f, Tuple{eltype(arg)})
    res = Vector{T}(n)
    nt = Threads.nthreads()
    Threads.@threads for _ in 1:nt
        while true
            i = Threads.atomic_add!(counter, 1)
            if i > n
                break
            end
            res[i] = f(arg[i])
        end
    end
    if !isleaftype(T)
        return [v for v in res]
    end
    return res
end

function create_sequence(t)
    builder = BuilderT(System.ThermalInit{1,Float32}(15, 4, 4), Setup.Dummy(),
                       System.HyperFineMeasure{3}())
    # builder = BuilderT(System.ThermalInit{1,Float32}(15, 4, 4), Setup.Dummy(),
    #                    Setup.CombinedMeasure(System.NBarMeasure(),
    #                                          System.GroundStateMeasure(),
    #                                          System.HyperFineMeasure{3}()))
    Setup.add_pulse(builder, System.RealRaman{Float32,1,3}(t, 1, (0f0, 0f0, 0f0),
                                                           (0, 0, 0), sz, []))
    return builder.seq
end

const xname = "t"
const params = linspace(0, 20, 100)
res = @time threadmap(p->Setup.run(create_sequence(p), statec(), nothing, 100000), params)

function plot_hook()
    xlabel(xname)
end

function plot_ground_state(params, res)
    figure()
    gp = [r.a for r in res]
    gp_unc = [r.s for r in res]
    errorbar(params, gp, gp_unc)
    title("Ground state probability")
    plot_hook()
    grid()
end

function plot_total(params, res)
    figure()
    total = [1 - t.a for t in res]
    total_unc = [t.s for t in res]
    errorbar(params, total, total_unc)
    title("Total loss")
    ylim([0, ylim()[2]])
    plot_hook()
    grid()
end

function plot_nbars(params, res)
    figure()
    nbarx_res = (nr[1] for nr in res)
    nbary_res = (nr[2] for nr in res)
    nbarz_res = (nr[3] for nr in res)

    nbarx = [n.a for n in nbarx_res]
    nbarx_unc = [n.s for n in nbarx_res]
    nbary = [n.a for n in nbary_res]
    nbary_unc = [n.s for n in nbary_res]
    nbarz = [n.a for n in nbarz_res]
    nbarz_unc = [n.s for n in nbarz_res]
    if max(maximum(nbarx), maximum(nbary), maximum(nbarz)) > 1
        gca()[:set_yscale]("log")
    end
    errorbar(params, nbarx, nbarx_unc, label="X")
    errorbar(params, nbary, nbary_unc, label="Y")
    errorbar(params, nbarz, nbarz_unc, label="Z")
    legend()
    title("\$\\bar n\$")
    plot_hook()
    grid()
end

plot_hf(params, res) = plot_hf(params, collect(res))

function plot_hf{T<:Tuple}(params, res::Vector{T})
    figure()
    for i in 1:nfields(T)
        hf = [r[i].a for r in res]
        hf_unc = [r[i].s for r in res]
        errorbar(params, hf, hf_unc, label="$i")
    end
    legend()
    ylim([0, 1])
    title("Hyperfine")
    plot_hook()
    grid()
end

function maybe_save(name)
    if !interactive()
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

function plot_result(params, res, prefix)
    plot_ground_state(params, (r[2] for r in res))
    maybe_save(joinpath(prefix, "ground.svg"))

    plot_total(params, (r[1][2] for r in res))
    maybe_save(joinpath(prefix, "loss.svg"))

    plot_nbars(params, (r[1][1] for r in res))
    maybe_save(joinpath(prefix, "nbar.svg"))

    plot_hf(params, (r[3][1] for r in res))
    maybe_save(joinpath(prefix, "hf.svg"))
end

# if interactive()
#     plot_result(params, res, "")
# else
#     mkpath(ARGS[1], 0o750)
#     plot_result(params, res, ARGS[1])
# end
plot_hf(params, [r[1] for r in res])
maybe_show()