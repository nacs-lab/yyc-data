#

__precompile__(true)

module NaCsPlot

import NaCsData
using NaCsCalc.Utils: interactive
using PyPlot
using PyCall

const rcParams = PyDict{PyAny,PyAny,true}(PyNULL())

function __init__()
    if !interactive()
        pygui(false)
    end
    copy!(PyObject(rcParams), matplotlib."rcParams")
    rcParams["svg.hashsalt"] = 19680801
    fontsize(20)
    ticksize(15)
    copy!(hist, PyPlot.matplotlib.pyplot.hist)
end

function nobold()
    rcParams["font.weight"] = "normal"
    return
end
function bold()
    rcParams["font.weight"] = "bold"
    return
end
function fontsize(s)
    rcParams["font.size"] = s
    return
end
function ticksize(s)
    matplotlib.rc("xtick", labelsize=s)
    matplotlib.rc("ytick", labelsize=s)
    return
end

const hist = PyCall.PyNULL()

function plot_data(data, columns, scale=1; yoffset=0, kws...)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    for col in columns
        errorbar(params, ratios[perm, col] .* scale .+ yoffset,
                 uncs[perm, col] .* scale; kws...)
    end
end

plot_loading_data(data, scale=1; yoffset=0, kws...) =
    plot_data(data, 1, scale; yoffset=yoffset, kws...)
plot_survival_data(data, scale=1; yoffset=0, kws...) =
    plot_data(data, 2, scale; yoffset=yoffset, kws...)

function save(name; close=true)
    dir = dirname(name)
    if !isempty(dir)
        mkpath(dir, mode=0o755)
    end
    metadata = Dict("Creator"=>nothing, "Producer"=>nothing, "CreationDate"=>nothing)
    savefig("$name.pdf"; bbox_inches="tight", transparent=true, metadata=metadata)
    savefig("$name.png"; bbox_inches="tight", transparent=true, metadata=metadata)
    savefig("$name.svg"; bbox_inches="tight", transparent=true,
            metadata=Dict("Creator"=>nothing, "Date"=>nothing,
                          "Format"=>nothing, "Type"=>nothing))
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
