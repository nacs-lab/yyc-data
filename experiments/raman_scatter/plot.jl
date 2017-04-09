#!/usr/bin/julia

include("load_data.jl")

using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const save_fig = get(ENV, "NACS_SAVE_FIG", "true") == "true"

function maybe_save(name)
    if save_fig
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if !save_fig
        show()
    end
end

const prefix = joinpath(@__DIR__, "imgs", "data_combined_20170405")

figure()
errorbar(params_f1_coprop, ratios_f1_coprop, uncs_f1_coprop, fmt="o-")
title("F1 co-prop")
grid()
ylim(0, ylim()[2])
xlim(0, xlim()[2])
xlabel("\$t (ms)\$")
maybe_save("$(prefix)_f1_coprop")

figure()
errorbar(params_f2_coprop, ratios_f2_coprop, uncs_f2_coprop, fmt="o-")
title("F2 co-prop")
grid()
ylim(0, ylim()[2])
xlim(0, xlim()[2])
xlabel("\$t (ms)\$")
maybe_save("$(prefix)_f2_coprop")

figure()
errorbar(params_f1_up, ratios_f1_up, uncs_f1_up, fmt="o-")
title("F1 up")
grid()
ylim(0, ylim()[2])
xlim(0, xlim()[2])
xlabel("\$t (ms)\$")
maybe_save("$(prefix)_f1_up")

figure()
errorbar(params_f1_down, ratios_f1_down, uncs_f1_down, fmt="o-")
title("F1 down")
grid()
ylim(0, ylim()[2])
xlim(0, xlim()[2])
xlabel("\$t (ms)\$")
maybe_save("$(prefix)_f1_down")

figure()
errorbar(params_f2_counterop, ratios_f2_counterop, uncs_f2_counterop, fmt="o-")
title("F2 counter-op")
grid()
ylim(0, ylim()[2])
xlim(0, xlim()[2])
xlabel("\$t (ms)\$")
maybe_save("$(prefix)_f2_counterop")

maybe_show()
