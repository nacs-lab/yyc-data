#!/usr/bin/julia -f

using PyPlot

data = readdlm("CsPhotonsVsDetuning-4.txt")
data2 = readdlm("CsPhotonsVsDetuning-7.txt")
data3 = readdlm("CsPhotonsVsDetuning-3.txt")

figure()
errorbar(data[:, 1], data[:, 2], data[:, 3],
         label="\$\\lambda=970nm\$", fmt="o-", linewidth=3, markersize=8)
errorbar(data2[:, 1], data2[:, 2], data2[:, 3],
         label="\$\\lambda=935nm\$", fmt="o-", linewidth=3, markersize=8)
ylim([0, 25000])
legend()
grid()
xlabel("Detuning (MHz)", size=20)
ylabel("Photon count", size=20)
# show()
savefig("det-vs-photon-dc.png", bbox_inches="tight")
close()

figure()
errorbar(data[:, 1], data[:, 2], data[:, 3],
         label="\$\\lambda=970nm\$", fmt="o-", linewidth=3, markersize=8)
errorbar(data2[:, 1], data2[:, 2], data2[:, 3],
         label="\$\\lambda=935nm\$", fmt="o-", linewidth=3, markersize=8)
ylim([0, 30000])
legend()
grid()
xlabel("Detuning (MHz)", size=20)
ylabel("Photon count", size=20)
# show()
savefig("det-vs-photon-dc2.png", bbox_inches="tight")
close()

figure()
errorbar(data[:, 1], data[:, 2], data[:, 3],
         label="\$\\lambda=970nm\$", fmt="o-", linewidth=3, markersize=8)
errorbar(data2[:, 1], data2[:, 2], data2[:, 3],
         label="\$\\lambda=935nm\$", fmt="o-", linewidth=3, markersize=8)
errorbar(data3[:, 1], data3[:, 2], data3[:, 3],
         label="\$\\lambda=970nm\$\n(modulated)", fmt="o-", linewidth=3, markersize=8)
ylim([0, 30000])
legend()
grid()
xlabel("Detuning (MHz)", size=20)
ylabel("Photon count", size=20)
# show()
savefig("det-vs-photon-ac.png", bbox_inches="tight")
close()
