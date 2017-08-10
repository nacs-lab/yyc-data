#!/usr/bin/julia

# Parameter, P_ground, σ(P_ground), loss, σ(loss),
# nbar_x, σ(nbar_x), nbar_y, σ(nbar_y), nbar_z, σ(nbar_z), hf1, σ(hf1), ..., hf8, σ(hf8)

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: interactive
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 25))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

function maybe_save(name)
    if !interactive()
        savefig("$name.pdf"; bbox_inches="tight", transparent=true)
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg"; bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

const datadir = joinpath(@__DIR__, "../../calculations/sideband/data")
const data_with_scatter = readcsv(joinpath(datadir, "qual_with_scatter.csv"), Float64)
const data_no_scatter = readcsv(joinpath(datadir, "qual_no_scatter.csv"), Float64)

figure()
errorbar(data_with_scatter[:, 1], data_with_scatter[:, 2], data_with_scatter[:, 3],
         label="With scatter", color="C1")
errorbar(data_no_scatter[:, 1], data_no_scatter[:, 2], data_no_scatter[:, 3],
         label="Without scatter", color="C0")
xlim([0, maximum(data_no_scatter[:, 1])])
ylim([0, 1])
xlabel("Cooling Cycle")
ylabel("Probability")
legend(fontsize=18, loc=4)
grid()
maybe_save(joinpath(@__DIR__, "imgs/simcool_real"))

figure()
errorbar(data_no_scatter[:, 1], data_no_scatter[:, 2], data_no_scatter[:, 3],
         label="Without scatter", color="C0")
xlim([0, maximum(data_no_scatter[:, 1])])
ylim([0, 1])
xlabel("Cooling Cycle")
ylabel("Probability")
legend(fontsize=18, loc=4)
grid()
maybe_save(joinpath(@__DIR__, "imgs/simcool_no_scatter"))

maybe_show()
