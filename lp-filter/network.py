#!/usr/bin/env python

from math import *
from impedance import Z, R, L, C
from numbers import Number

# u_out(i_out) --z1--> ----> u_in(i_in)
#                     |
#                     z2
#                     |
def _transform_z12(u, i, z1, z2, n=1):
    for _ in range(n):
        i = i + u / z2
        u =  u + z1 * i
    return u, i

class Network(object):
    def __init__(self, node, n=1):
        self.__node = node
        self.__n = n
    def z(self, f, load):
        for i in range(self.__n):
            load = self.__node.z(f, load)
        return load
    def transform(self, f, u_in, i_in, n=1):
        return self.__node.transform(f, u_in, i_in, n * self.__n)
    def gain(self, f, load):
        if load is None:
            u = 1
            i = 0
        else:
            if not isinstance(load, Number):
                load = load.z(f)
            u = 1
            i = 1 / load
        u, i = self.transform(f, u, i)
        return 1 / abs(u)

class NetworkNode(Network):
    def __init__(self, z1, z2):
        if isinstance(z1, Number):
            z1 = R(z1)
        if isinstance(z2, Number):
            z2 = R(z2)
        self.__z1 = z1
        self.__z2 = z2
    def z(self, f, load):
        if load is None:
            return (self.__z1 + self.__z2).z(f)
        elif load == 0:
            return self.__z1.z(f)
        return (self.__z1 + (self.__z2 & load)).z(f)
    def transform(self, f, u_in, i_in, n=1):
        return _transform_z12(u_in, i_in, self.__z1.z(f), self.__z2.z(f), n)

class RCLPFilter(Network):
    def __init__(self, r, c, n=1):
        Network.__init__(self, NetworkNode(R(r), C(c)), n)

class RCHPFilter(Network):
    def __init__(self, r, c, n=1):
        Network.__init__(self, NetworkNode(C(c), R(r)), n)

class LCLPFilter(Network):
    def __init__(self, l, c, n=1):
        Network.__init__(self, NetworkNode(L(l), C(c)), n)

class LCHPFilter(Network):
    def __init__(self, l, c, n=1):
        Network.__init__(self, NetworkNode(C(c), L(l)), n)
