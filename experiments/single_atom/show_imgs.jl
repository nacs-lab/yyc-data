#!/usr/bin/julia

using MAT
using PyPlot

# PyPlot init
PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 30,
                                            "font.weight" => "bold"))
PyPlot.matplotlib[:rc]("xtick", labelsize=25)
PyPlot.matplotlib[:rc]("ytick", labelsize=25)

# Load image
imgs = matopen(ARGS[1]) do fd
    read(fd, "Images")
end
imgs = imgs[:, :, 1:2:end]

# Single atom
roix, roiy = 19:22, 19:22
roi_img = imgs[roix, roiy, :]
counts = sum(roi_img, [1, 2])[1, 1, :]
single_atom = counts .>= 38

# Plot single atom
# figure()
# PyPlot.matplotlib[:pyplot][:hist](counts, bins=40)

# Substract background
atom_imgs = imgs[:, :, single_atom]
no_atom_imgs = imgs[:, :, (!).(single_atom)]
no_bgd = atom_imgs .- mean(no_atom_imgs, 3)[:, :, 1]

# Plot smaller sub image
# sub_img = imgs[20:21, 20:21, :]
# figure()
# PyPlot.matplotlib[:pyplot][:hist](sum(sub_img, [1, 2])[1, 1, :], bins=40)

# Plot various average image
# figure()
# implot = imshow(mean(atom_imgs, 3)[11:31, 11:31, 1], interpolation="none")
# implot[:set_cmap]("viridis")

# figure()
# implot = imshow(mean(no_bgd, 3)[11:31, 11:31, 1], interpolation="none")
# implot[:set_cmap]("viridis")

# figure()
# implot = imshow(mean(atom_imgs, 3)[19:22, 19:22, 1], interpolation="none")
# implot[:set_cmap]("viridis")

std_img = std(no_bgd, 3)[:, :, 1]
avg_img = mean(no_bgd, 3)[:, :, 1]
figure()
imshow(std_img[11:31, 11:31], interpolation="none", cmap="viridis")
colorbar()
figure()
imshow(avg_img[11:31, 11:31], interpolation="none", cmap="viridis")
colorbar()

function calc_corrs(imgs, xy, regions)
    x, y = xy
    roix, roiy = regions
    # For each pair of values, calculate (<v1 v2> - <v1> <v2>) / (N)
    nx = length(roix)
    ny = length(roiy)
    s12 = zeros(eltype(imgs), nx, ny)
    s2 = zeros(eltype(imgs), nx, ny)
    s1 = 0
    nimg = size(imgs, 3)
    @inbounds for idx in 1:nimg
        v1 = imgs[x, y, idx]
        s1 += v1
        for iy in 1:ny
            ry = roiy[iy]
            for ix in 1:nx
                rx = roix[ix]
                v2 = imgs[rx, ry, idx]
                s2[ix, iy] += v2
                s12[ix, iy] += v2 * v1
            end
        end
    end
    a1 = s1 / nimg
    # @inbounds for i in 1:length(s12)
    #     s12[i] = s12[i] / s2[i] / a1 - 1
    # end
    @inbounds for i in 1:length(s12)
        s12[i] = s12[i] - a1 * s2[i]
    end
    scale!(s12, 1 / maximum(s12))
    return s12
end
# corrs20_20 = calc_corrs(imgs, (20, 20), (11:31, 11:31))
# figure()
# imshow(corrs20_20, interpolation="none", cmap="bwr", vmin=-1, vmax=1)
# colorbar()

# corrs24_24 = calc_corrs(imgs, (24, 24), (11:31, 11:31))
# figure()
# imshow(corrs24_24, interpolation="none", cmap="bwr", vmin=-1, vmax=1)
# colorbar()

show()
