using Scical
using PyPlot
using LsqFit

function tryParseNumber(s)
    try
        parse(Float64, s)
    catch
        s
    end
end

function getFileInfo(name)
    bname = basename(name)
    return ("$name.txt", Any[tryParseNumber(v) for v in split(bname, '_')])
end

function getData(finfo, unc, offset)
    # 15% + 0.02μW uncertainty by default
    data = readdlm(finfo[1], ' ', Float64)
    idxs = sortperm(data[:, 1])
    res = Array{Float64, 2}(size(data, 1), 3)
    res[:, 1:2] = data[idxs, :]
    for i in 1:size(data, 1)
        res[i, 3] = res[i, 2] * unc + offset
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

function plotSingleData(finfo, data)
    figure()
    plotData(finfo, data)
    xlim([minimum(data[:, 1]) - 5, maximum(data[:, 1]) + 5])
    ylim([0, maximum(data[:, 2]) * 1.1])
end

model(x, p) = begin
    angl = x / 180 * π
    angl2 = 2 * angl
    p[1] + p[2] * sin(angl2) + p[3] * cos(angl2)
end

function fitData(finfo, data)
    fit = curve_fit(model, data[:, 1], data[:, 2], 1 ./ data[:, 3].^2,
                    [mean(data[:, 2]), 0.0, 0.0])

    # We can estimate errors on the fit parameters,
    # to get 95% confidence error bars:
    errors = estimate_errors(fit, 0.95)
    (fit.param, errors)
end

function plotDataFit(finfo, data)
    param, errors = fitData(finfo, data)
    plotSingleData(finfo, data)
    xs = linspace(minimum(data[:, 1]) - 5, maximum(data[:, 1]) + 5, 1000)
    ys = model(xs, param)
    plot(xs, ys)
end

plotFileFit(finfo, unc=0.01, offset=0.02e-3) =
    plotDataFit(finfo, getData(finfo, unc, offset))
