#!/usr/bin/julia

using PyPlot
PyPlot.matplotlib["rcParams"][:update](Dict("font.size" => 20))
PyPlot.matplotlib[:rc]("xtick", labelsize=15)
PyPlot.matplotlib[:rc]("ytick", labelsize=15)

X_orders = [1, -1, -2]
X_freqs = [18.23, 19.196, 19.650]
X_uncs = [0.010, 0.010, 0.010]

Y_orders = [1, -1, -2]
Y_freqs = [18.22, 19.221, 19.714]
Y_uncs = [0.010, 0.010, 0.010]

Z_orders = [1, -1, -2, -3, -4, -5]
Z_freqs = [18.5325, 18.708, 18.793, 18.884, 18.973, 19.055]
Z_uncs = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]

errorbar(X_orders, X_freqs, X_uncs, label="X")
errorbar(Y_orders, Y_freqs, Y_uncs, label="Y")
errorbar(Z_orders, Z_freqs, Z_uncs, label="Z")

show()
