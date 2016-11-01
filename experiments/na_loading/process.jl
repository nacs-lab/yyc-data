#!/usr/bin/julia

using PyPlot
using NaCsData

matplotlib["rcParams"][:update](Dict("font.size" => 30,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=25)
matplotlib[:rc]("ytick", labelsize=25)

const params, ratios, uncs = NaCsData.calc_survival(ARGS[1])

const freqs = linspace(-5e6, -30e6, 11)
const base_amps = [0, 0.685, 1.32, 1.82, 2.36, 2.75, 3.24,
                   3.75, 4.29, 4.83, 5.30, 5.84, 6.41]
const amp_scals = [6.62, 6.56, 6.40, 6.32, 6.21, 5.98,
                   5.70, 5.40, 5.01, 4.46, 3.97] / 6.40
const nfreqs = length(freqs)
const namps = length(base_amps)
const all_freqs = Vector{Float64}(nfreqs * namps)
const all_amps = Vector{Float64}(nfreqs * namps)

for i in 1:nfreqs
    local i
    local idx_offset = (i - 1) * namps
    local freq = freqs[i]
    local amp_scal = amp_scals[i]
    for j in 1:namps
        local j
        idx = j + idx_offset
        all_freqs[idx] = freq
        all_amps[idx] = base_amps[j] * amp_scal
    end
end

const ratios_matrix = reshape(ratios[:, 1], namps, nfreqs)

# The following is a hacky 2D interpolation implementation since I can't
# get any non-hacky one to work quickly.....
# It split the 2D surface into triangles and do linear interpolation in each
# triangle.
function dist_to_line(p1, p2, xy)
    x1, y1 = p1
    x2, y2 = p2
    x0, y0 = xy
    y21 = y2 - y1
    x21 = x2 - x1
    return abs(y21 * x0 - x21 * y0 + x2 * y1 - y2 * x1) / sqrt(y21^2 + x21^2)
end

function calc_val(f, a)
    fi = (f - first(freqs)) / step(freqs)
    fi < 0 && return NaN
    fi_lo = floor(Int, fi) + 1
    fi_hi = fi_lo + 1
    fi_hi > nfreqs && return NaN
    fx = fi - fi_lo + 1
    amp_scal_lo = amp_scals[fi_lo]
    amp_scal_hi = amp_scals[fi_hi]

    # Figure out which square this fall into
    aj_lo = 0
    aj_hi = 0
    right_f = false
    for j in 2:namps
        base_a = base_amps[j]
        a_lo = base_a * amp_scal_lo
        a_hi = base_a * amp_scal_hi
        if a_lo * (1 - fx) + a_hi * fx >= a
            aj_lo = j - 1
            aj_hi = j
            right_f = base_amps[j - 1] * amp_scal_lo * (1 - fx) + a_hi * fx >= a
            break
        end
    end
    aj_lo == 0 && return NaN
    xy = (a, fx)
    if right_f
        p1 = base_amps[aj_lo] * amp_scal_hi, 1.0
        p2 = base_amps[aj_hi] * amp_scal_hi, 1.0
        p3 = base_amps[aj_lo] * amp_scal_lo, 0.0
        v1 = ratios_matrix[aj_lo, fi_hi]
        v2 = ratios_matrix[aj_hi, fi_hi]
        v3 = ratios_matrix[aj_lo, fi_lo]
    else
        p1 = base_amps[aj_lo] * amp_scal_lo, 0.0
        p2 = base_amps[aj_hi] * amp_scal_lo, 0.0
        p3 = base_amps[aj_hi] * amp_scal_hi, 1.0
        v1 = ratios_matrix[aj_lo, fi_lo]
        v2 = ratios_matrix[aj_hi, fi_lo]
        v3 = ratios_matrix[aj_hi, fi_hi]
    end
    w1 = dist_to_line(p2, p3, xy) / dist_to_line(p2, p3, p1)
    w2 = dist_to_line(p1, p3, xy) / dist_to_line(p1, p3, p2)
    w3 = dist_to_line(p2, p1, xy) / dist_to_line(p2, p1, p3)
    return w1 * v1 + w2 * v2 + w3 * v3
end

const freqs_p = linspace(minimum(freqs), maximum(freqs), 2000)
const amps_p = linspace(0, maximum(amp_scals) * maximum(base_amps), 2000)

plot_img = [calc_val(f, a) for a in amps_p, f in freqs_p]

figure()
imshow(plot_img, vmin=0.5,
       extent=[freqs_p[1] / 1e6, freqs_p[end] / 1e6, amps_p[end], amps_p[1]],
       aspect="auto")
xlabel("Detuning (MHz)")
ylabel("Power (mW)")
colorbar()
# savefig(ARGS[2], bbox_inches="tight", transparent=true)
show()
