#!/usr/bin/julia

include("process.jl")

const prefix = ARGS[1]

# Name, Power, File
const num_lines = (length(ARGS) - 1) ÷ 3
@assert length(ARGS) == (num_lines * 3 + 1)

figure()
for i in 1:num_lines
    name = ARGS[i * 3 - 1]
    powers = readcsv(ARGS[i * 3])
    params, ratios, uncs = NaCsData.calc_survival(ARGS[i * 3 + 1])

    ((depths, ratios1, uncs1),
     (plot_depth, fit_cdf, fit_pdf)) = process_lowering(powers, params,
                                                        ratios, uncs)
    plot(plot_depth, fit_pdf, label=name)
end

ylim([0, ylim()[2]])
xlabel("\$T\$ (MHz)")
ylabel("Population")
legend()
grid()
# show()
savefig("$(prefix).png", bbox_inches="tight", transparent=true)
close()
