using Scical
using PyPlot
using LsqFit

function getFileInfo(name)
    bname = basename(name)
    return ("$name.txt", Any[parse(Float64, v) for v in split(bname, '_')])
end

function getData(finfo)
    data = readdlm(finfo[1], ' ', Float64)
    idxs = sortperm(data[:, 1])
    res = Array{Float64, 2}(size(data, 1), 3)
    res[:, 1:2] = data[idxs, :]
    for i in 1:size(data, 1)
        res[i, 3] = res[i, 2] * 0.01 + 0.02e-3 # 15% + 0.02μW
    end
    return res
end
function plotData(finfo, data)
    errorbar(data[:, 1], data[:, 2], data[:, 3])
    lamb, waveplate = finfo[2][1:2]
    nplate = if waveplate == 0
        "None"
    elseif waveplate == 4
        "Quarter waveplate"
    elseif waveplate == 2
        "Half waveplate"
    else
        error("Unknown waveplate.")
    end
    angl_str = if waveplate != 0
        angl = finfo[2][3]
        " @ \$$angl^\\circ\$"
    else
        ""
    end
    title("\$$(Int(lamb))\$nm, $nplate$angl_str")
end

function plotFile(finfo)
    data = getData(finfo)
    plotData(finfo, data)
    data
end

function plotSingleData(finfo, data)
    figure()
    plotData(finfo, data)
    xlim([minimum(data[:, 1]) - 5, maximum(data[:, 1]) + 5])
    ylim([0, maximum(data[:, 2]) * 1.1])
end

function plotSingleFile(finfo)
    plotSingleData(finfo, getData(finfo))
end
