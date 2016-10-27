#!/usr/bin/julia

using MAT
using PyPlot

PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 30,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

imgs = matopen(ARGS[1]) do fd
    read(fd, "Images")
end
imgs = imgs[:, :, 1:2:end]

roi = 19:22, 19:22
roi_img = imgs[roi..., :]
figure()
count = sum(roi_img, [1, 2])[1, 1, :]
single_atom = count .>= 38
PyPlot.matplotlib[:pyplot][:hist](count, bins=40)

atom_imgs = imgs[:, :, single_atom]
no_atom_imgs = imgs[:, :, (!).(single_atom)]
no_bgd = max.(atom_imgs .- mean(no_atom_imgs, 3)[:, :, 1], 0)

roi_img = imgs[20:21, 20:21, :]
figure()
count = sum(roi_img, [1, 2])[1, 1, :]
single_atom = count .>= 38
PyPlot.matplotlib[:pyplot][:hist](count, bins=40)

figure()
implot = imshow(mean(atom_imgs, 3)[11:31, 11:31, 1], interpolation="none")
implot[:set_cmap]("viridis")

figure()
implot = imshow(mean(no_bgd, 3)[11:31, 11:31, 1], interpolation="none")
implot[:set_cmap]("viridis")

figure()
implot = imshow(mean(atom_imgs, 3)[19:22, 19:22, 1], interpolation="none")
implot[:set_cmap]("viridis")
show()
