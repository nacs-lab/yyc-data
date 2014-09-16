#!/usr/bin/env python

from numpy import *

def sequence_gen(mod_freq, mod_amp, cent_freq, dt=100e-6, T=1):
    t = -dt
    while t < T:
        t += dt
        yield cent_freq + mod_amp * sin(mod_freq * 2 * pi * t)

def main():
    print('dt = 10 us, TTL(all) = 0')
    for freq in sequence_gen(2, 20e6, 80e6, dt=1e-3, T=.5):
        print('dt = 1000 us, freq(2) = %.0f Hz' % freq)

if __name__ == '__main__':
    main()
