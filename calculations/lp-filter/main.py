#!/usr/bin/env python

from pylab import *
from network import Network, NetworkNode, RCLPFilter
from impedance import Z, R, L, C

def calc_gains(fs, r, c, n):
    rclp = RCLPFilter(r, c, n)
    return [rclp.gain(f, None) for f in fs]

def main():
    fs = r_[0:2e4:2000j][1:]
    plot(log10(fs), log10(calc_gains(fs, 2, 1e-6, 4)) * 10, label="4")
    grid()
    legend()
    show()

if __name__ == '__main__':
    main()
