#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsPlot
using PyPlot

include("clock_shift.jl")

const prefix = joinpath(@__DIR__, "imgs/real0")

const rg0 = load_dir(joinpath(@__DIR__, "data/real0/0"))
const rg1 = load_dir(joinpath(@__DIR__, "data/real0/2"))
δs = get_δ(rg0)
states = 1:9

const plt_data_dir = joinpath(@__DIR__, "plot_data")
mkpath(plt_data_dir, 0o755)

function write_dlmstr(io, dlm, vals; print_first=false)
    first = !print_first
    for v in vals
        if !first
            print(io, dlm, v)
        else
            print(io, v)
            first = false
        end
    end
end

write_energycsv(fname, rg) = open(fname, "w") do io
    write(io, "State")
    write_dlmstr(io, ',', get_δ(rg), print_first=true)
    println(io)
    for i in 1:50
        write_dlmstr(io, ' ', rg.h.states[i])
        write_dlmstr(io, ',', get_energy(rg, i), print_first=true)
        println(io)
    end
end

write_overlapcsv(fname, rg, init) = open(fname, "w") do io
    write(io, "State")
    write_dlmstr(io, ',', get_δ(rg), print_first=true)
    println(io)
    for i in 1:50
        write_dlmstr(io, ' ', rg.h.states[i])
        write_dlmstr(io, ',', get_overlap(rg, i, init), print_first=true)
        println(io)
    end
end

figure()
for i in states
    plot(δs ./ 1000, get_energy(rg0, i) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Energy (kHz)")
ylim([-100, 100])
title("Parity $(rg0.p) 9 lowest")
grid()
NaCsPlot.maybe_save("$(prefix)_energy_000_9")

figure()
for i in states
    plot(δs ./ 1000, get_energy(rg1, i) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Energy (kHz)")
ylim([50, 200])
title("Parity $(rg1.p) 9 lowest")
grid()
NaCsPlot.maybe_save("$(prefix)_energy_010_9")

init_state0 = get_state(rg0, (0, 0, 0, 0, 0, 0))
states0 = collect(filter_overlap(rg0, init_state0, 0.3))
figure()
for i in states0
    plot(δs ./ 1000, get_energy(rg0, i) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Energy (kHz)")
ylim([-100, 155])
title("Parity $(rg0.p) ⩾ 0.3 overlap")
grid()
NaCsPlot.maybe_save("$(prefix)_energy_000_0")

figure()
for i in states0
    plot(δs ./ 1000, get_overlap(rg0, i, init_state0))
end
xlabel("1st order shift (kHz)")
ylabel("Overlap")
title("Parity $(rg0.p) ⩾ 0.3 overlap")
ylim([0, 1])
grid()
NaCsPlot.maybe_save("$(prefix)_overlap_000_0")

init_state1 = get_state(rg1, (0, 0, 0, 0, 1, 0))
states1 = collect(filter_overlap(rg1, init_state1, 0.3))
figure()
for i in states1
    plot(δs ./ 1000, get_energy(rg1, i) ./ 1000)
end
xlabel("1st order shift (kHz)")
ylabel("Energy (kHz)")
ylim([50, 200])
title("Parity $(rg1.p) ⩾ 0.3 overlap")
grid()
NaCsPlot.maybe_save("$(prefix)_energy_010_0")

figure()
for i in states1
    plot(δs ./ 1000, get_overlap(rg1, i, init_state1))
end
xlabel("1st order shift (kHz)")
ylabel("Overlap")
ylim([0, 1])
title("Parity $(rg1.p) ⩾ 0.3 overlap")
grid()
NaCsPlot.maybe_save("$(prefix)_overlap_010_0")

NaCsPlot.maybe_show()

write_energycsv(joinpath(plt_data_dir, "energy_000.csv"), rg0)
write_energycsv(joinpath(plt_data_dir, "energy_010.csv"), rg1)
write_overlapcsv(joinpath(plt_data_dir, "overlap_000.csv"), rg0, init_state0)
write_overlapcsv(joinpath(plt_data_dir, "overlap_010.csv"), rg1, init_state1)
