push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
import NaCsCalc: Trap
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)
