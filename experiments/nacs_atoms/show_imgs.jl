#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using MAT
using PyPlot
using Statistics

# NaCsPlot.fontsize(30)
# NaCsPlot.ticksize(25)
# NaCsPlot.bold()

const iname = joinpath(@__DIR__, "data", "data_20161019_171310.mat")

imgs, single_atom = matopen(iname) do fd
    read(fd, "Images"), read(fd, "SingleAtom") .!= 0
end

const whitelist = [27, 55, 78, 93, 144, 154, 171, 245, 311, 329, 542, 622, 629, 631, 632, 678,
                   897, 917, 918, 933, 961, 962, 1067, 1101, 1160, 1237, 1393, 1430, 1465,
                   1468, 1486, 1490, 1645, 1695, 1729, 1719]

const prefix = joinpath(@__DIR__, "imgs")

figure()
implot = imshow(imgs[11:31, 11:31, 1695], interpolation="nearest",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("gray")
xlabel("x(\$\\mathrm{\\mu m}\$)")
ylabel("y(\$\\mathrm{\\mu m}\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
NaCsPlot.maybe_save("$(prefix)/single_gray")

figure()
implot = imshow(imgs[11:31, 11:31, 1695], interpolation="nearest",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("viridis")
xlabel("x(\$\\mathrm{\\mu m}\$)")
ylabel("y(\$\\mathrm{\\mu m}\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
NaCsPlot.maybe_save("$(prefix)/single_viridis")

figure()
implot = imshow(mean(imgs, dims=3)[11:31, 11:31, 1], interpolation="nearest",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("gray")
xlabel("x(\$\\mathrm{\\mu m}\$)")
ylabel("y(\$\\mathrm{\\mu m}\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
NaCsPlot.maybe_save("$(prefix)/avg_gray")

figure()
implot = imshow(mean(imgs, dims=3)[11:31, 11:31, 1], interpolation="nearest",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("viridis")
xlabel("x(\$\\mathrm{\\mu m}\$)")
ylabel("y(\$\\mathrm{\\mu m}\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
NaCsPlot.maybe_save("$(prefix)/avg_viridis")

function show_imgs(imgs, single_atom)
    i = 1
    nimgs = size(imgs, 3)
    while i <= nimgs
        fig = figure()
        for j in 1:25
            while i <= nimgs
                i in whitelist && break
                # single_atom[i] && break
                i += 1
            end
            i <= nimgs || break
            fig[:add_subplot](5, 5, j)
            implot = imshow(imgs[10:31, 10:31, i], interpolation="nearest")
            implot[:set_cmap]("gray")
            # implot[:set_cmap]("viridis")
            title(i)
            i += 1
        end
    end
end

if interactive() && false
    show_imgs(imgs, single_atom)
end

NaCsPlot.maybe_show()
