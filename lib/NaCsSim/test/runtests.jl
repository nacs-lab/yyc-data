#!/usr/bin/julia

@everywhere module TestSequence
import NaCsSim: Setup, System
import NaCsCalc: Trap

const BuilderT = Setup.SeqBuilder{System.StateC,Void}
const sz = 500, 100, 100

const m_Na = 23f-3 / 6.02f23
const k_Na = Float32(2π) / 589f-9
η(freq) = Trap.η(m_Na, freq, k_Na)

const η_op = (η(60f3), η(420f3), η(420f3))
const η_op_dri = (0f0, η(420f3) / sqrt(2f0), η(420f3) / sqrt(2f0))
const isσ = [false false false
              true true true
              true true true]
const branching_11 = Float32[0 0 1 / 6
                             0 0 1 / 12
                             0 0 1 / 4]
const branching_21 = Float32[0 1 / 3 0
                             0 1 / 6 0
                             0 1 / 2 0]
const branching_22 = Float32[1 / 3 0 0
                             1 / 6 0 0
                             1 / 2 0 0]
const η_raman1 = (η(60f3) / sqrt(2f0), 0f0, 0f0)
const η_raman2 = (0f0, η(420f3) * sqrt(2f0), 0f0)
const η_raman3 = (0f0, 0f0, η(420f3) * sqrt(2f0))
const η_ramans = (η_raman1, η_raman2, η_raman3)

# 1: (2, -2); 2: (2, -1); 3: (1, -1)
export statec
const statec = System.StateC(sz...)
function op_pulse(t, γ1, γ2, γ2′=0.01f0 * γ2, ηs=η_op, ηdri=η_op_dri)
    branching = (branching_11 .* γ1 .+ branching_21 .* γ2 .+
                 branching_22 .* γ2′)
    System.OP{Float32}(t, branching, ηs, ηdri, isσ)
end
function raman_pulse(ax, order, t, Ω=1)
    ns = ntuple(i->i == ax ? -order : 0, Val{3})
    System.Raman{Float32,1,3}(t, Ω, η_ramans[ax], ns, sz)
end

@inline function mix(p, delta, i)
    return p + delta * (i - 1)
end

immutable OPParams
    γ1::Float32
    γ2::Float32
    darkness::Float32
end
immutable OPDelta
    γ1::Float32
    γ2::Float32
    darkness::Float32
end
Base.:*(delta::OPDelta, n) = OPDelta(delta.γ1 * n, delta.γ2 * n, delta.darkness * n)
Base.:+(param::OPParams, delta::OPDelta) =
    OPParams(param.γ1 + delta.γ1, param.γ2 + delta.γ2,
             param.darkness + delta.darkness)
OPDelta(;γ1=0, γ2=0, darkness=0) = OPDelta(γ1, γ2, darkness)

pulse(params::OPParams) =
    op_pulse(1, params.γ1, params.γ2, params.γ2 * params.darkness)

immutable RamanParams
    ax::Int
    order::Int
    t::Float32
end
immutable RamanDelta
    t::Float32
end
Base.:*(delta::RamanDelta, n) = RamanDelta(delta.t * n)
Base.:+(param::RamanParams, delta::RamanDelta) =
    RamanParams(param.ax, param.order, param.t + delta.t)

pulse(params::RamanParams) = raman_pulse(params.ax, params.order, params.t)

add_pulse(builder, params) = Setup.add_pulse(builder, pulse(params))

# Group with two axial cooling per loop
immutable Grp2AParams
    op::OPParams
    raman11::RamanParams
    raman12::RamanParams
    raman2::RamanParams
    raman3::RamanParams
    deltaop::OPDelta
    delta11::RamanDelta
    delta12::RamanDelta
    delta2::RamanDelta
    delta3::RamanDelta
    rep::Int
end

function add_raman_op(builder, raman, op)
    raman.t > 0 || return
    add_pulse(builder, raman)
    add_pulse(builder, op)
    return
end

function add_pulse(builder, params::Grp2AParams)
    n = params.rep
    if n == 0
        return
    end
    for i in 1:n
        op = params.op + params.deltaop * i
        raman11 = params.raman11 + params.delta11 * i
        raman12 = params.raman12 + params.delta12 * i
        raman2 = params.raman2 + params.delta2 * i
        raman3 = params.raman3 + params.delta3 * i
        add_raman_op(builder, raman11, op)
        add_raman_op(builder, raman12, op)
        add_raman_op(builder, raman2, op)
        add_raman_op(builder, raman11, op)
        add_raman_op(builder, raman12, op)
        add_raman_op(builder, raman3, op)
    end
end

# TODO sweep OP power
export create_sequence
function create_sequence(ncycles)
    # ncycles = 88
    op_defect = 0.01
    pulses_left = Ref(ncycles)
    cooling_on = true

    take_pulses = function (n)::Int
        cooling_on || return 0
        nleft = pulses_left[]
        if nleft > n
            nleft -= n
            pulses_left[] = nleft
            return n
        else
            pulses_left[] = 0
            return nleft
        end
    end

    builder = BuilderT(System.ThermalInit{1,Float32}(15, 4, 4), Setup.Dummy(),
                       Setup.CombinedMeasure(System.NBarMeasure(),
                                             System.GroundStateMeasure(),
                                             System.HyperFineMeasure{3}()))

    pulses = [Grp2AParams(OPParams(15, 0.8, op_defect),
                          RamanParams(1, 6, 5),
                          RamanParams(1, 5, 5),
                          RamanParams(2, 2, 4),
                          RamanParams(3, 2, 4),
                          OPDelta(),
                          RamanDelta(0 / 11),
                          RamanDelta(5 / 11),
                          RamanDelta(2 / 11),
                          RamanDelta(2 / 11),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.7, op_defect),
                          RamanParams(1, 5, 5),
                          RamanParams(1, 4, 5),
                          RamanParams(2, 2, 4),
                          RamanParams(3, 2, 4),
                          OPDelta(),
                          RamanDelta(5 / 11),
                          RamanDelta(5 / 11),
                          RamanDelta(2 / 11),
                          RamanDelta(2 / 11),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.7, op_defect),
                          RamanParams(1, 4, 5),
                          RamanParams(1, 3, 4.5),
                          RamanParams(2, 2, 5),
                          RamanParams(3, 2, 5),
                          OPDelta(),
                          RamanDelta(5 / 11),
                          RamanDelta(4.5 / 11),
                          RamanDelta(7 / 11),
                          RamanDelta(7 / 11),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.4, op_defect),
                          RamanParams(1, 3, 10),
                          RamanParams(1, 2, 5),
                          RamanParams(2, 1, 3),
                          RamanParams(3, 1, 3),
                          OPDelta(γ2=0.2 / 11),
                          RamanDelta(0 / 11),
                          RamanDelta(5 / 11),
                          RamanDelta(2 / 11),
                          RamanDelta(2 / 11),
                          take_pulses(12)),
              Grp2AParams(OPParams(15, 0.3, op_defect),
                          RamanParams(1, 2, 5),
                          RamanParams(1, 1, 1.5),
                          RamanParams(2, 1, 2.5),
                          RamanParams(3, 1, 2.5),
                          OPDelta(γ2=-0.08 / 11),
                          RamanDelta(0 / 39),
                          RamanDelta(1 / 39),
                          RamanDelta(2 / 39),
                          RamanDelta(2 / 39),
                          take_pulses(40))
              ]
    for p in pulses
        add_pulse(builder, p)
    end
    return builder.seq
end

end

@everywhere import NaCsSim: Setup, System
@everywhere using TestSequence

# const params = linspace(0.0, 0.02, 41)
const params = 0:88
const xname = "Cycles"

res = pmap(p->Setup.run(create_sequence(p), statec, nothing, 100000), params)

using PyPlot
PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 15,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

function plot_ground_state(params, res)
    figure()
    gp = [r.a for r in res]
    gp_unc = [r.s for r in res]
    errorbar(params, gp, gp_unc)
    title("Ground state probability")
    xlabel(xname)
    grid()
end

function plot_total(params, res)
    figure()
    total = [1 - t.a for t in res]
    total_unc = [t.s for t in res]
    errorbar(params, total, total_unc)
    title("Total loss")
    xlabel(xname)
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
    xlabel(xname)
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
    title("Hyperfine")
    xlabel(xname)
    grid()
end

function plot_result(params, res)
    plot_ground_state(params, (r[2] for r in res))
    savefig(joinpath(ARGS[1], "ground.svg"), bbox_inches="tight", transparent=true)
    savefig(joinpath(ARGS[1], "ground.png"), bbox_inches="tight", transparent=true)

    plot_total(params, (r[1][2] for r in res))
    savefig(joinpath(ARGS[1], "loss.svg"), bbox_inches="tight", transparent=true)
    savefig(joinpath(ARGS[1], "loss.png"), bbox_inches="tight", transparent=true)

    plot_nbars(params, (r[1][1] for r in res))
    savefig(joinpath(ARGS[1], "nbar.svg"), bbox_inches="tight", transparent=true)
    savefig(joinpath(ARGS[1], "nbar.png"), bbox_inches="tight", transparent=true)

    plot_hf(params, (r[3][1] for r in res))
    savefig(joinpath(ARGS[1], "hf.svg"), bbox_inches="tight", transparent=true)
    savefig(joinpath(ARGS[1], "hf.png"), bbox_inches="tight", transparent=true)
end

plot_result(params, res)
show()
