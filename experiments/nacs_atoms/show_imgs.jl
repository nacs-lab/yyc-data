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

# implot = imshow(imgs[11:31, 11:31, 1695], interpolation="none",
#                 extent=[-2.5, 2.5, -2.5, 2.5])
# implot[:set_cmap]("gray")
# # implot[:set_cmap]("viridis")
# xlabel("x(\$\\mu m\$)")
# ylabel("y(\$\\mu m\$)")
# savefig(ARGS[2], bbox_inches="tight")
# show()

implot = imshow(mean(imgs, 3)[11:31, 11:31, 1], interpolation="none",
                extent=[-2.5, 2.5, -2.5, 2.5])
# implot[:set_cmap]("gray")
implot[:set_cmap]("viridis")
xlabel("x(\$\\mu m\$)")
ylabel("y(\$\\mu m\$)")
savefig(ARGS[2], bbox_inches="tight")
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
