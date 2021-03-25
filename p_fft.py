########################################################################
# program: fft.py
# author: Tom Irvine
# Email: tom@vibrationdata.com
# version: 1.9
# date: September 12, 2013
# description:
#
#  Determine the Fast Fourier transform of a signal.
#  The file must have two columns: time(sec) & amplitude.
#
########################################################################

from __future__ import print_function

import sys
from tkinter import Label

# if sys.version_info[0] == 2:
#     print ("Python 2.x")
#     import Tkinter as tk
#     from tkFileDialog import asksaveasfilename


if sys.version_info[0] == 3:
    print("Python 3.x")
    import tkinter as tk
    from tkinter.filedialog import asksaveasfilename


from p_tompy import read_two_columns_from_dialog, signal_stats, sample_rate_check
from p_tompy import GetInteger2, WriteData3
from p_tompy import time_history_plot

from math import pi, atan2, log
from numpy import zeros, argmax, mean

from scipy.fftpack import fft

import matplotlib.pyplot as plt


########################################################################

class FFT:

    def __init__(self, a, b, imr):
        self.a = a
        self.b = b
        self.imr = imr
        pass

    def fft_data(self):

        #   Truncate to 2**n

        num = len(self.b)

        noct = int(log(num)/log(2))

        num_fft = 2**noct

        bb = self.b[0:num_fft]

        if(self.imr == 1):
            bb = bb-mean(bb)

        dur_fft = self.a[num_fft-1]-self.a[0]

        df = 1/dur_fft

        z = fft(bb)

        nhalf = num_fft//2

        print(" ")
        print(" %d samples used for FFT " % num_fft)
        print("df = %8.4g Hz" % df)

        zz = zeros(nhalf, 'f')
        ff = zeros(nhalf, 'f')
        ph = zeros(nhalf, 'f')

        freq = zeros(num_fft, 'f')

        z /= float(num_fft)

        for k in range(0, int(num_fft)):
            freq[k] = k*df

        ff = freq[0:nhalf]

        for k in range(0, int(nhalf)):

            if(k > 0):
                zz[k] = 2.*abs(z[k])
            else:
                zz[k] = abs(z[k])

            ph[k] = atan2(z.real[k], z.imag[k])

        idx = argmax(abs(zz))

        return idx, freq, ff, z, zz, ph, nhalf, df, num_fft


########################################################################

class READ_DATA_FFT:

    def __init__(self):
        pass

    @classmethod
    def check_data(cls, a, b, num, sr, dt):

        sample_rate_check(a, b, num, sr, dt)

        return sr, dt

    def read_and_stats(self):

        # label = "Enter the time history filename"
        label = "acceleration_data.csv"

        a, b, num = read_two_columns_from_dialog(label)

        sr, dt, ave, sd, rms, skew, kurtosis, dur = signal_stats(a, b, num)

        sr, dt = READ_DATA_FFT.check_data(a, b, num, sr, dt)

        return a, b, num, sr, dt, dur, ave, rms

#######################################################################
