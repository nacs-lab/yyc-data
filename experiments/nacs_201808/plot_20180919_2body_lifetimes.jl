#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

load_mat(name) = NaCsData.load_striped_mat(joinpath(@__DIR__, "data", name))

names = [
    "data_20180919_230048.mat",
    "data_20180920_000758.mat",
    "data_20180920_011532.mat",
    "data_20180920_022920.mat",
    "data_20180920_033832.mat",
    "data_20180920_044858.mat",
    "data_20180920_055809.mat"
]
files = [load_mat(name) for name in names]

split_data(spec, f, ls...) =
    NaCsData.split_data(NaCsData.select_count(f..., NaCsData.select_single(ls...)), spec)

const spec_a = OrderedDict(
    :cool=>[0, 10, 20, 50, 100, 200, 400],
    :nocool=>[0, 10, 20, 50, 100, 200, 400],
)

freqs = 690:696

datas_na = [split_data(spec_a, f, (1, -2), (3, -4)) for f in files]
datas_cs = [split_data(spec_a, f, (-1, 2), (-3, 4)) for f in files]
datas_nacs_na = [split_data(spec_a, f, (1, 2), (3,)) for f in files]
datas_nacs_cs = [split_data(spec_a, f, (1, 2), (4,)) for f in files]

const prefix = joinpath(@__DIR__, "imgs", "data_20180919_2body_lifetimes")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_na[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_cs[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Cs, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_na[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_nocool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_cs[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Cs, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_nocool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_na[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na + Cs, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_na_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_cs[i][:cool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Na + Cs, cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_cs_cool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_na[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Na / Na + Cs, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_na_nocool")

figure()
for i in 1:7
    NaCsPlot.plot_survival_data(datas_nacs_cs[i][:nocool], fmt=".-", label="$(freqs[i])")
end
grid()
legend()
xlim([0, 500])
ylim([0, 1])
title("Cs / Na + Cs, not cooled")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_cs_nocool")

# lifetime_model(x, p) = p[1] .* exp.(.-x ./ p[2])
# lifetime_model2(x, p) = p[1] .* exp.(.-x ./ p[2]) .+ p[3]

# const plotx = linspace(0, 0.5, 1001)
# fit_na = fit_survival(lifetime_model, data_na_lifetime, [0.95, 1.0], plotx=plotx)
# fit_cs = fit_survival(lifetime_model, data_cs_lifetime, [0.95, 1.0], plotx=plotx)
# fit_nacs = fit_survival(lifetime_model, data_nacs_lifetime, [0.95, 0.2], plotx=plotx)
# fit_nacs2 = fit_survival(lifetime_model2, data_nacs_lifetime, [0.95, 0.2, 0.1], plotx=plotx)

# function fit_text2(xy1, xy2, name, fit; kws...)
#     t = Unc(fit.param[2], fit.unc[2])
#     offset = Unc(fit.param[3], fit.unc[3])
#     text(xy1..., "\$\\tau_{$name}=$(t)s\$"; kws...)
#     text(xy2..., "\$off=$(offset)\$"; kws...)
# end

# figure()
# NaCsPlot.plot_survival_data(data_cs_lifetime, fmt="C0.")
# plot(fit_cs.plotx, fit_cs.ploty, "C0-", label="Cs")
# NaCsPlot.plot_survival_data(data_na_lifetime, fmt="C1.")
# plot(fit_na.plotx, fit_na.ploty, "C1-", label="Na")
# NaCsPlot.plot_survival_data(data_nacs_lifetime, fmt="C2.")
# plot(fit_nacs.plotx, fit_nacs.ploty, "C2-", label="Na+Cs")
# plot(fit_nacs2.plotx, fit_nacs2.ploty, "C3-", label="With offset")
# text(0.3, 0.9, "\$\\tau_{Cs}=$(Unc(fit_cs.param[2], fit_cs.unc[2]))s\$", color="C0")
# text(0.05, 0.7, "\$\\tau_{Na}=$(Unc(fit_na.param[2], fit_na.unc[2]))s\$", color="C1")
# text(0.03, 0.05, "\$\\tau_{Na+Cs}=$(Unc(fit_nacs.param[2], fit_nacs.unc[2]))s\$", color="C2")
# fit_text2((0.32, 0.25), (0.32, 0.17), "", fit_nacs2, color="C3")
# grid()
# legend(fontsize=16)
# xlim([0, 0.55])
# ylim([0, 1])
# title("Lifetime")
# xlabel("Time (s)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
