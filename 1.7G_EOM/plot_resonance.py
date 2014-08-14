#!/usr/bin/env python

import matplotlib
matplotlib.rc('font', size=20, weight='bold')

from scipy.optimize import curve_fit, leastsq
from pylab import *
import inspect

class Ret(dict):
    def __init__(self, *args, **kwargs):
        for arg in args:
            try:
                self._a(**arg)
            except:
                frame = inspect.currentframe().f_back
                try:
                    self[arg] = frame.f_locals[arg]
                except:
                    self[arg] = frame.f_globals[arg]
        self._a(**kwargs)
    def _a(self, **kwargs):
        self.update(kwargs)
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)
    def __setattr__(self, key, value):
        self[key] = value
    def __delattr__(self, name):
        del self[name]
    def __getitem_iter__(self, keys):
        for k in keys:
            yield dict.__getitem__(self, k)
    def __getitem__(self, keys):
        if isscalar(keys):
            return dict.__getitem__(self, keys)
        return self.__getitem_iter__(keys)
    def __setitem__(self, keys, items):
        if isscalar(keys):
            dict.__setitem__(self, keys, items)
            return
        for k, v in zip(keys, items):
            dict.__setitem__(self, k, v)
    def __dir__(self):
        return list(self.keys()) + dict.__dir__(self)

def redchi2(delta, sigma, n):
    '''chi2 / dof'''
    return sum((delta / sigma)**2) / (delta.size - n)

def curve_fit_wrapper(fitfun):
    # http://www.physics.utoronto.ca/~phy326/python/curve_fit_to_data.py
    # Notes: maxfev is the maximum number of func evaluations tried; you
    # can try increasing this value if the fit fails.
    # If the program returns a good chi-squared but an infinite
    # covariance and no parameter uncertainties, it may be
    # because you have a redundant parameter;
    # try fitting with a simpler function.
    def curve_fitter(x, y, sig=None, p0=None):
        a, cov = curve_fit(fitfun, x, y, sigma=sig, p0=p0)
        yfit = fitfun(x, *a)
        func = lambda x: fitfun(x, *a)
        return Ret('x', 'a', 'yfit', 'cov', 'func', s=sqrt(diag(cov)),
                   chi2=redchi2(y - yfit, sig, len(a)) if sig != None else None)
    return curve_fitter

def a_pm_s(a_s, unit='', sci=None, tex=False):
    try:
        a = a_s.a
        s = a_s.s
    except AttributeError:
        try:
            a = a_s[0]
            s = a_s[1]
        except KeyError:
            a = a_s['a']
            s = a_s['s']

    try:
        a = [i for i in a]
        l = len(a)
    except TypeError:
        a = [a]
        l = 1

    try:
        s = [i for i in s]
    except TypeError:
        s = [s]

    if len(s) < l:
        s += [0] * (l - len(s))

    if type(unit) == type(''):
        unit = [unit] * l
    else:
        try:
            unit = [u for u in unit]
        except TypeError:
            unit = [unit] * l

    if len(unit) < l:
        unit += [''] * (l - len(unit))
    if l == 1:
        return _a_pm_s(a[0], s[0], unit[0], sci, tex)
    return array([_a_pm_s(a[i], s[i], unit[i], sci, tex) for i in range(0, l)])

def _a_pm_s(a, s, unit, sci, tex):
    '''input: observable,error
       output: formatted observable +- error in scientific notation'''
    if s <= 0:
        return '%f%s' % (a, unit)

    if sci == None:
        if s < 100 and (abs(a) > 1 or s > 1):
            sci = False
        else:
            sci = True

    la = int(floor(log10(abs(a))))
    ls = int(floor(log10(s)))
    fs = floor(s * 10**(1 - ls))
    if sci:
        fa = a * 10**-la
        dl = la - ls + 1
    else:
        fa = a
        dl = 1 - ls
    dl = dl if dl > 0 else 0

    if dl == 1:
        ss = '%.1f' % (fs / 10)
    else:
        ss = '%.0f' % fs

    if sci:
        if tex:
            return (('%.' + ('%d' % dl) + r'f(%s)\times10^{%d}{%s}') %
                    (fa, ss, la, unit))
        else:
            return ('%.' + ('%d' % dl) + 'f(%s)*10^%d%s') % (fa, ss, la, unit)
    else:
        return ('%.' + ('%d' % dl) + 'f(%s)%s') % (fa, ss, unit)

freqs, main1, main2, side1, side2 = loadtxt('resonance-2014_8_14.txt').T
sigma_y = 3
offset_y = 4

data_l = len(freqs)

main1 -= offset_y
main2 -= offset_y
side1 -= offset_y
side2 -= offset_y

main1_s = ones(data_l) * sigma_y
main2_s = ones(data_l) * sigma_y
side1_s = ones(data_l) * sigma_y
side2_s = ones(data_l) * sigma_y

ratio1 = side1 / main1 * 100
ratio2 = side2 / main2 * 100
ratio1_s = sqrt((side1_s / side1)**2 + (main1_s / main1)**2) * ratio1
ratio2_s = sqrt((side2_s / side2)**2 + (main2_s / main2)**2) * ratio2

@curve_fit_wrapper
def fit_resonance(x, x0, y0, w):
    return y0 / (1 + (2 * (x - x0) / w)**2)

freqs_all = r_[freqs, freqs]
ratio_all = r_[ratio1, ratio2]
ratio_all_s = r_[ratio1_s, ratio2_s]
res = fit_resonance(freqs_all, ratio_all, ratio_all_s,
                    p0=[mean(freqs), max(ratio_all), std(freqs)])

freq_min = min(freqs)
freq_max = max(freqs)

freqs_smooth = r_[freq_min - 1:freq_max + 1:1000j]
yfit_smooth = res.func(freqs_smooth)

# print('x0', a_pm_s([res.a[0], res.s[0]], 'MHz', tex=True))
# print('y0', a_pm_s([res.a[1], res.s[1]], '\\%', tex=True))
# print('w', a_pm_s([res.a[2], res.s[2]], 'MHz', tex=True))

errorbar(freqs, ratio1, ratio1_s, fmt='bo')
errorbar(freqs, ratio2, ratio2_s, fmt='bo')
plot(freqs_smooth, yfit_smooth, 'r')
xlim([min(freqs) - 1, max(freqs) + 1])
xlabel('Frequency (MHz)')
ylabel('Ratio of sideband to main peak (%)')
title('$1.7\\mathbf{GHz}$ EOM resonance')
text(min(freqs) - .8, ylim()[1],
     '$R = \\frac{R_{max}}{1 + \\left(2\\cdot(f - f_0) / \Gamma\\right)^2}$\n'
     '$f_0 = %s$\n'
     '$R_{max} = %s$\n'
     '$\Gamma = %s$\n' % (a_pm_s([res.a[0], res.s[0]],
                                 '\\mathbf{MHz}', tex=True),
                          a_pm_s([res.a[1], res.s[1]],
                                 '\\mathbf{\\%}', tex=True),
                          a_pm_s([res.a[2], res.s[2]],
                                 '\\mathbf{MHz}', tex=True)),
     horizontalalignment='left', verticalalignment='top')

grid()
savefig('resonance-2014_8_14.png')
# show()
