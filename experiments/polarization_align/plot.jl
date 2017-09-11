#!/usr/bin/julia

using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 10,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

errorbar([-4, -6, -8, -10, -12, -14] - 12 * 2,
         [0.830, 0.730, 0.623, 0.640, 0.620, 0.796],
         [0.056, 0.033, 0.039, 0.039, 0.042, 0.032], label="H=12\$^\\circ\$")
errorbar([-4, -6.5, -8, -10] - 14 * 2,
         [0.623, 0.578, 0.603, 0.758],
         [0.042, 0.046, 0.027, 0.053], label="H=14\$^\\circ\$")
errorbar([-10, -8, -6, -4, -2, 0, 2] - 16 * 2,
         [0.816, 0.806, 0.691, 0.573, 0.557, 0.657, 0.713],
         [0.040, 0.037, 0.033, 0.039, 0.041, 0.038, 0.037], label="H=16\$^\\circ\$")
errorbar([4, 2, 0, -2] - 18 * 2,
         [0.656, 0.554, 0.567, 0.664],
         [0.037, 0.036, 0.040, 0.037], label="H=18\$^\\circ\$")
errorbar([0, 2, 4, 6] - 20 * 2,
         [0.745, 0.563, 0.562, 0.664],
         [0.043, 0.040, 0.041, 0.032], label="H=20\$^\\circ\$")
errorbar([12, 10, 8] - 22 * 2,
         [0.772, 0.807, 0.736],
         [0.035, 0.038, 0.035], label="H=22\$^\\circ\$")

grid()
legend()
xlabel("Q - 2H (\$^\\circ\$)", fontsize=20)
ylabel("F1 population", fontsize=20)
ylim([0, 0.9])
savefig("2017-09-11.pdf", bbox_inches="tight", transparent=true)
savefig("2017-09-11.png", bbox_inches="tight", transparent=true)
