#!/usr/bin/env python

from __future__ import division

from numpy import *

freq0 = 1070.5971e6
d_freq = 30e6
ramp_len = 1000e-6

def gen_freq(start, end, t_len):
    num = int(t_len / 10e-6)
    diff = end - start
    return ((start + diff * i / (num - 1)) for i in range(num))

def main():
    print('dt = 10 us, TTL(all) = 0')
    print('dt = 10 us, TTL(16) = 1')
    for freq in gen_freq(freq0, freq0 + d_freq, ramp_len):
        print('dt = 10 us, freq(19) = %.0f Hz' % freq)
    print('dt = 10 us, TTL(17) = 1')
    print('dt = 10000 us, TTL(16) = 0')
    print('dt = 10 us, TTL(17) = 0')
    for freq in gen_freq(freq0 + d_freq, freq0, ramp_len):
        print('dt = 10 us, freq(19) = %.0f Hz' % freq)
    print('dt = 10000 us, TTL(all) = 0')

if __name__ == '__main__':
    main()
