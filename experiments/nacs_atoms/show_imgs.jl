#!/usr/bin/julia

using MAT
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 30))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

imgs, single_atom = matopen(ARGS[1]) do fd
    read(fd, "Images"), read(fd, "SingleAtom") .!= 0
end

const whitelist = [27, 55, 78, 93, 144, 154, 171, 245, 311, 329, 542, 622, 629, 631, 632, 678, 897, 917, 918, 933, 961, 962, 1067, 1101, 1160, 1237, 1393, 1430, 1465, 1468, 1486, 1490, 1645, 1695, 1729, 1719]

implot = imshow(imgs[11:31, 11:31, 1695], interpolation="none",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("gray")
xlabel("x(\$\\mu m\$)")
ylabel("y(\$\\mu m\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
savefig(joinpath(ARGS[2], "single_gray.png"),
        bbox_inches="tight", transparent=true)
close()

implot = imshow(imgs[11:31, 11:31, 1695], interpolation="none",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("viridis")
xlabel("x(\$\\mu m\$)")
ylabel("y(\$\\mu m\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
savefig(joinpath(ARGS[2], "single_viridis.png"),
        bbox_inches="tight", transparent=true)
close()

implot = imshow(mean(imgs, 3)[11:31, 11:31, 1], interpolation="none",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("gray")
xlabel("x(\$\\mu m\$)")
ylabel("y(\$\\mu m\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
savefig(joinpath(ARGS[2], "avg_gray.png"),
        bbox_inches="tight", transparent=true)
close()

implot = imshow(mean(imgs, 3)[11:31, 11:31, 1], interpolation="none",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("viridis")
xlabel("x(\$\\mu m\$)")
ylabel("y(\$\\mu m\$)")
annotate("Na", xy=(-0.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2, 0.5), xytext=(2.5, -4),
         arrowprops=Dict("facecolor"=>"white", "shrink"=>0.1,
                         "width"=>15, "headwidth"=>25, "headlength"=>25),
         color="white", size=40, weight="bold")
savefig(joinpath(ARGS[2], "avg_viridis.png"),
        bbox_inches="tight", transparent=true)
close()
# show()

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
    show()
end

# show_imgs(imgs, single_atom)
