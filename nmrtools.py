#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nmrtools

Tools for analyzing long runs of pulse NMR data to find the Larmor
frequency and T2* for each pulse. The incoming data are expected to
be 2-D arrays where arr[0, :] are the times of the samples and arr[1, :]
the actual samples. The expectation is that the data will consist of a
sequence of pulses separated by periods of silence.

Created on Thu Oct 19 11:08:33 2023

@author: bcollett
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def decsin(t, a, f, t0, phase=0):
    # print('ds: ', a, f, t0, phase)
    s = np.sin(2 * np.pi * f * (t - t[0]) + phase)
    return a * np.exp(-(t - t[0])/t0) * s


#
# A Pulse is just a short stream of data, broken out into separate time
# and value arrays. It knows how to usine fitting to extract the best-fit
# parameters. It also knows how to plot itself, with or without the fit.
#
# You need a set of starting parameters for a successful curve_fit on
# complex data like this. I make the starting parameters class variables
# so that if you want to alter them then you can just re-assign them.
#
class Pulse:
    # Class Variables
    init_amp = 0.7
    init_freq = 2000
    init_t2 = 0.0042
    init_phase = 2.5

    def __init__(self, t, v):
        self.times = t
        self.vals = v
        self.amp = self.init_amp
        self.freq = self.init_freq
        self.t2 = self.init_t2
        self.phase = self.init_phase

    def plot(self, addFit=False):
        plt.plot(self.times, self.vals, '.', markersize=2)
        if addFit:
            print(self.amp, self.freq, self.t2, self.phase)
            fv = decsin(self.times, self.amp, self.freq, self.t2, self.phase)
            plt.plot(self.times, fv)

    def fit(self, initp=None):
        if initp is None:
            initp = [self.init_amp, self.init_freq,
                     self.init_t2, self.init_phase]
            print(initp)
        popt, pcov = curve_fit(decsin, self.times, self.vals, p0=initp)
        # print(popt)
        self.amp = popt[0]
        self.freq = popt[1]
        self.t2 = popt[2]
        self.phase = popt[3]
        return popt, pcov


#
# I use a PulseIterator to split the stream of data into the pulses it
# contains, discarding the silence between. Note that the each pulse
# retains its full time so that all the times for pulse n are greater
# than all the times for pulse n-1.
#
# The pulse detection process is inherently a little complex and so I
# have made it depend on a set of parameters that are class variables.
# If you want non-default values then simply re-assign the class variables
# before creating the iterator.
#
class PulseIterator:
    # Class Variables
    max_silence = 0.2
    pulse_length = 0.025

    def __init__(self, stream):
        self.times = stream[0, :]
        self.values = stream[1, :]
        plt.figure()
        plt.title('Data', loc='center')
        plt.plot(self.times, self.values)
        self._start = 0
        self._end = len(self.times)
        self._npp = int(self.pulse_length/(self.times[1]-self.times[0]))

    def __iter__(self):
        return self

    def __next__(self):
        #
        # Run in loop skipping values < max and aborting if fall off end.
        # print(self._end)
        while True:
            # print(self._start)
            if self._start >= self._end:
                raise StopIteration
            if np.abs(self.values[self._start]) > self.max_silence:
                # Found start of a pulse
                end = self._start + self._npp
                if end > self._end:
                    end = self._end
                print(self._start, self._end)
                pulse_times = self.times[self._start:end]
                pulse_vals = self.values[self._start:end]
                self._start = end
                return Pulse(pulse_times, pulse_vals)
            else:
                self._start += 1

def extract(stream):
    pit = PulseIterator(stream)
    pulses = []
    for p in pit:
        pulses.append(p)
    return pulses

def getFreqs(pulses, ip=[0.4, 2000, 0.002, 1.8]):
    pulse_t = []
    pulse_f = []
    print(ip)
    for p in pulses:
        fp, fc = p.fit(ip)
        print(p.times[0], fp[1])
        pulse_t.append(p.times[0])
        pulse_f.append(fp[1])
    return pulse_t, pulse_f

if __name__ == '__main__':
    strm = np.loadtxt('Fids100s.txt')
#    pt, pv, pits = extract(strm)
#    plt.plot(pt, pv, '.')
    ps = extract(strm)
    plt.figure()
    plt.title('Pulses', loc='center')
    for p in ps:
        plt.plot(p.times, p.vals, '.')
    pts, pfs = getFreqs(ps)
    plt.figure()
    plt.plot(pts, pfs, '+')