#!/usr/bin/env python

from math import *
from numbers import Number

class Impedance(object):
    def z(self, f):
        raise NotImplementedError
    def __add__(self, other):
        if not isinstance(other, Impedance) or isinstance(other, Number):
            return NotImplemented
        return Series(self, other)
    def __radd__(self, other):
        if not isinstance(other, Impedance) or isinstance(other, Number):
            return NotImplemented
        return Series(other, self)
    def __and__(self, other):
        if not isinstance(other, Impedance) or isinstance(other, Number):
            return NotImplemented
        return Parrallel(self, other)
    def __rand__(self, other):
        if not isinstance(other, Impedance) or isinstance(other, Number):
            return NotImplemented
        return Parrallel(other, self)

Z = Impedance

class Series(Z):
    def __init__(self, z1, z2):
        if isinstance(z1, Number):
            z1 = R(z1)
        if isinstance(z2, Number):
            z2 = R(z2)
        self.__z1 = z1
        self.__z2 = z2
    def z(self, f):
        return self.__z1.z(f) + self.__z2.z(f)

class Parallel(Z):
    def __init__(self, z1, z2):
        if isinstance(z1, Number):
            z1 = R(z1)
        if isinstance(z2, Number):
            z2 = R(z2)
        self.__z1 = z1
        self.__z2 = z2
    def z(self, f):
        return 1 / (1 / self.__z1.z(f) + 1 / self.__z2.z(f))

class Resistor(Z):
    def __init__(self, r):
        self.__r = r
    def z(self, f):
        return self.__r

R = Resistor

class Capacitor(Z):
    def __init__(self, c):
        self.__c = c
    def z(self, f):
        return 1 / (1j * 2 * pi * f * self.__c)

C = Capacitor

class Inductor(Z):
    def __init__(self, l):
        self.__l = l
    def z(self, f):
        return 1j * 2 * pi * f * self.__l

L = Inductor
