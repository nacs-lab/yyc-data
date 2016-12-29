#!/usr/bin/julia

using MAT
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 30,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)
const pyhist = PyPlot.matplotlib[:pyplot][:hist]

imgs, single_atom = matopen(ARGS[1]) do fd
    read(fd, "Images"), read(fd, "SingleAtom") .!= 0
end

const na_roi = 20:22, 20:21
const cs_roi = 19:20, 25:26

function hist_for_roi(imgs, roi)
    sub_imgs = imgs[roi..., :]
    counts = sum(sub_imgs, (1, 2))[1, 1, :]
    pyhist(counts, bins=30)
end

figure()
hist_for_roi(imgs[:, :, 2:2:end], cs_roi)
title("Cesium")
xlabel("Electron count")
savefig(joinpath(ARGS[2], "hist_cs.svg"),
        bbox_inches="tight", transparent=true)
close()

figure()
hist_for_roi(imgs[:, :, 2:2:end], na_roi)
title("Sodium")
xlabel("Electron count")
savefig(joinpath(ARGS[2], "hist_na.svg"),
        bbox_inches="tight", transparent=true)
close()

# show()
