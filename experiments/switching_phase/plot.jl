#!/usr/bin/julia

using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

data = sortrows(readcsv(ARGS[1], header=true)[1])
phases = data[:, 1]
loading = data[:, 2]
loading_unc = sqrt.(loading .* 10) ./ 10
avg1 = data[:, 3]
avg1_unc = data[:, 4]
avg2 = data[:, 5]
avg2_unc = data[:, 6]

total_count = avg1 .- avg2
total_count_unc = sqrt.(avg1_unc.^2 + avg2_unc.^2)
avg_count = total_count ./ loading .* 100
avg_count_unc = sqrt.((loading_unc ./ loading).^2 + (total_count_unc ./ total_count).^2) .* avg_count

fig, ax1 = subplots()
p1 = ax1[:errorbar](phases, loading, loading_unc, label="Loading", fmt="r.-")
ax1[:set_xlabel]("Phase (ns)")
ax1[:set_ylabel]("Loading (%)", color="r")
ax1[:tick_params]("y", colors="r")

ax2 = ax1[:twinx]()
p2 = ax2[:errorbar](phases, total_count .* 2, total_count_unc .* 2, label="Total count", fmt="b.-")
p3 = ax2[:errorbar](phases, avg_count, avg_count_unc, label="Average count", fmt="g.-")
ax2[:set_ylabel]("Count", color="b")
ax2[:tick_params]("y", colors="b")

ax1[:legend]([p1, p2, p3], ["Loading", "Total count", "Average count"], loc=4,
             fontsize=15)
fig[:tight_layout]()
grid()
savefig("$(ARGS[2])/combined.svg", bbox_inches="tight", transparent=true)
savefig("$(ARGS[2])/combined.png", bbox_inches="tight", transparent=true)
# show()
