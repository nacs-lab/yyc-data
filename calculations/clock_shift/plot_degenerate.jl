#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsPlot
using PyPlot

include("clock_shift.jl")

const prefix = joinpath(@__DIR__, "imgs/degenerate")

const rg0 = load_dir(joinpath(@__DIR__, "data/degenerate/0"))
const rg1 = load_dir(joinpath(@__DIR__, "data/degenerate/4"))
δs = get_δ(rg0)
states = 1:16

figure()
for i in states
    plot(δs ./ 1000, get_energy(rg0, i) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Energy (kHz)")
title("Parity $(rg0.p)")
grid()
NaCsPlot.maybe_save("$(prefix)_energy_000")

figure()
for i in states
    plot(δs ./ 1000, get_energy(rg1, i) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Energy (kHz)")
title("Parity $(rg1.p)")
grid()
NaCsPlot.maybe_save("$(prefix)_energy_100")

figure()
for i in states
    plot(δs ./ 1000, get_energy(rg0, i, base=rg0.r0.vals[i]) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Shift (kHz)")
title("Parity $(rg0.p)")
grid()
NaCsPlot.maybe_save("$(prefix)_shift_000")

figure()
for i in states
    plot(δs ./ 1000, get_energy(rg1, i, base=rg1.r0.vals[i]) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Shift (kHz)")
title("Parity $(rg1.p)")
grid()
NaCsPlot.maybe_save("$(prefix)_shift_100")

NaCsPlot.maybe_show()
