#

__precompile__(true)

module NaCsPlot

import NaCsData
using NaCsCalc.Utils: interactive
using PyPlot

function __init__()
    matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                         "font.weight" => "bold"))
    matplotlib[:rc]("xtick", labelsize=15)
    matplotlib[:rc]("ytick", labelsize=15)
end

function plot_survival_data(data, scale=1; kws...)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = ratios[perm, 2] .* scale
    uncs = uncs[perm, 2] .* scale
    errorbar(params, ratios, uncs; kws...)
end

function save(name; close=true)
    dir = dirname(name)
    if !isempty(dir)
        mkpath(dir, 0o755)
    end
    savefig("$name.pdf"; bbox_inches="tight", transparent=true)
    savefig("$name.png"; bbox_inches="tight", transparent=true)
    savefig("$name.svg"; bbox_inches="tight", transparent=true)
    close && PyPlot.close()
    return
end

function maybe_save(name)
    if !interactive()
        save(name)
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

end
