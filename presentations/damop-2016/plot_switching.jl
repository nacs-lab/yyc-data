#!/usr/bin/julia -f

using PyPlot

ts = linspace(0, 4π, 10000)
trap = atan(exp(15 * cos(ts) - 10)) * 2
img = atan(exp(-20 * cos(ts))) * 2

figure()
plot(ts, trap, "b", label="Trap", linewidth=3)
plot(ts, img, "r", label="Resonant", linewidth=3)
ylim([0, 4])
axis("off")
legend()
# show()
savefig("switching.png", bbox_inches="tight")
