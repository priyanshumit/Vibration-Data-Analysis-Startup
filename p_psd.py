########################################################################
# program: psd.py
# author: Tom Irvine
# version: 1.7
# date: September 12, 2013
# description:
#
#  Determine the power spectral density of a signal.
#  The file must have two columns: time(sec) & amplitude.
#
########################################################################

from __future__ import print_function
import datetime
import os
from scipy.signal import find_peaks
import numpy as np

import sys

# if sys.version_info[0] == 2:
#     print("Python 2.x")
#     import Tkinter as tk
#     from tkFileDialog import asksaveasfilename


if sys.version_info[0] == 3:
    print("Python 3.x")
    import tkinter as tk
    from tkinter.filedialog import asksaveasfilename


from p_tompy import read_two_columns_from_dialog, signal_stats, sample_rate_check
from p_tompy import GetInteger2, WriteData2
from p_tompy import time_history_plot

from sys import stdin
from math import sqrt, pi, log
from numpy import zeros, argmax, linspace, cos, mean
from scipy.fftpack import fft

import matplotlib.pyplot as plt

########################################################################


class READ_DATA_PSD:

    def __init__(self):
        pass

    @classmethod
    def check_data(cls, a, b, num, sr, dt):

        sample_rate_check(a, b, num, sr, dt)

        return sr, dt

    def read_and_stats(self):

        # label = "Enter the input time history"
        label = "acceleration_data.csv"

        a, b, num = read_two_columns_from_dialog(label)

        sr, dt, ave, sd, rms, skew, kurtosis, dur = signal_stats(a, b, num)

        sr, dt = READ_DATA_PSD.check_data(a, b, num, sr, dt)

        return a, b, num, sr, dt, dur

########################################################################


def GetString():
    iflag = 0
    while iflag == 0:
        try:
            s = stdin.readline()
            iflag = 1
        except ValueError:
            print('Invalid String')
    return s

########################################################################


def select_options(num, dt, window):

    # print(" ")
    # print(" Remove mean:  1=yes  2=no ")

    # mr_choice = GetInteger2()
    mr_choice = 1

    # print(" ")
    # print(" Select Window: 1=Rectangular 2=Hanning ")

    if window == "Rectangular":
        h_choice = 1
    elif window == "Hann":
        h_choice = 2
    elif window == "Blackman":
        h_choice = 3
    elif window == "Hamming":
        h_choice = 4

    # h_choice = GetInteger2()

    n = num

    ss = zeros(n)
    seg = zeros(n, 'f')
    i_seg = zeros(n)

    NC = 0
    for i in range(0, 1000):
        nmp = 2**(i-1)
        if(nmp <= n):
            ss[i] = 2**(i-1)
            seg[i] = float(n)/float(ss[i])
            i_seg[i] = int(seg[i])
            NC = NC+1
        else:
            break

    print(' ')
    print(' Number of   Samples per   Time per        df    ')
    print(' Segments     Segment      Segment(sec)   (Hz)   dof')

    for i in range(1, NC+1):
        j = NC+1-i
        if j > 0:
            if(i_seg[j] > 0):
                tseg = dt*ss[j]
                ddf = 1./tseg
                print('%8d \t %8d \t %10.3g  %10.3g    %d'
                      % (i_seg[j], ss[j], tseg, ddf, 2*i_seg[j]))
        if(i == 12):
            break

    ijk = 0
    while ijk == 0:
        # print(' ')
        # print(' Choose the Number of Segments:  ')
        # s = stdin.readline()
        s = 1
        NW = int(s)
        for j in range(0, len(i_seg)):
            if NW == i_seg[j]:
                ijk = 1
                break

# check

    mmm = 2**int(log(float(n)/float(NW))/log(2))

    df = 1./(mmm*dt)

# begin overlap

    mH = ((mmm/2)-1)

    return mmm, NW, df, mH, mr_choice, h_choice


########################################################################

class PSD:

    def __init__(self, mmm, NW, mH, df, b, mr_choice, h_choice):
        self.mmm = mmm
        self.NW = NW
        self.mH = mH
        self.df = df
        self.b = b
        self.mr_choice = mr_choice
        self.h_choice = h_choice

    @classmethod
    def Hanning_initial(cls, mmm):
        H = zeros(mmm, 'f')
        tpi = 2*pi
        alpha = linspace(0, tpi, mmm)
        ae = sqrt(8./3.)
        H = ae*0.5*(1.-cos(alpha))
        return H

    @classmethod
    def magnitude_resolve(cls, mmm, mH, Y):
        #
        mHm1 = mH-1
        z = zeros(int(mH), 'f')
        mag_seg = zeros(int(mH), 'f')
#
        z = abs(Y)/float(mmm)
#
        mag_seg[0] = z[0]**2
#
        mag_seg[1:int(mHm1)] = ((2*z[1:int(mHm1)])**2)/2
#
        return mag_seg

    def psd_core(self):

        if self.h_choice == 2:
            H = PSD.Hanning_initial(self.mmm)

        print(" ")
        print("     number of segments   NW= %d " % self.NW)
        print("       samples/segments  mmm= %d " % self.mmm)
        print(" half samples/segment-1   mH=%d  " % self.mH)
        print(" ")
        print("        df=%6.3f Hz" % self.df)

        full = zeros(int(self.mH), 'f')
        mag_seg = zeros(int(self.mH), 'f')

        amp_seg = zeros(self.mmm, 'f')

        nov = 0

        for ijk in range(1, 2*self.NW):

            amp_seg[0:self.mmm] = self.b[(0+nov):(self.mmm+nov)]

            nov = nov+int(self.mmm/2)

            if self.mr_choice == 1 or self.h_choice == 2:
                amp_seg -= mean(amp_seg)

            if self.h_choice == 2:
                amp_seg *= H

            Y = fft(amp_seg)

            mag_seg = PSD.magnitude_resolve(self.mmm, self.mH, Y)

            full += mag_seg

        den = self.df*(2*self.NW-1)

        full /= den

        ms = sum(full)

        freq = zeros(int(self.mH), 'f')

        maxf = (self.mH-1)*self.df

        freq = linspace(0, maxf, int(self.mH))

        tempf = freq[0:int(self.mH)-1]
        tempa = full[0:int(self.mH)-1]
        freq = tempf
        full = tempa

        rms = sqrt(ms*self.df)

        return rms, freq, full

########################################################################


def psd_plots(a, b, freq, full, rms, idx):

    pmin = 10**40
    pmax = 10**-40

    fmin = 10**40
    fmax = 10**-40

    for i in range(0, len(freq)):
        if full[i] > 0 and freq[i] > 0 and full[i] > pmax:
            pmax = full[i]
        if full[i] > 0 and freq[i] > 0 and full[i] < pmin:
            pmin = full[i]
        if freq[i] > 0 and freq[i] > fmax:
            fmax = freq[i]
        if freq[i] > 0 and freq[i] < fmin:
            fmin = freq[i]

    xmax = 10**-30
    xmin = xmax

    for i in range(-30, 30):
        if(fmax < 10**i):
            xmax = 10**i
            break

    for i in range(30, -30, -1):
        if(fmin > 10**i):
            xmin = 10**i
            break

    ymax = 10**-30
    ymin = ymax

    for i in range(-30, 30):
        if(pmax < 10**i):
            ymax = 10**i
            break

    for i in range(30, -30, -1):
        if(pmin > 10**i):
            ymin = 10**i
            break

    # print(" ")
    # print(" Is the input data dimension Accel(G) ?")
    # print(" 1=yes  2=no")

    # ind = GetInteger2()
    ind = 1

    if(ind == 1):
        th_label = 'Accel (G)'
        psd_label = 'Accel (G^2/Hz)'
    else:
        print('Enter input amplitude unit ')
        th_label = GetString()
        print('Enter PSD unit label, i.e.  unit^2/Hz')
        psd_label = GetString()

    print(" ")
    print(" view plots ")

    time_history_plot(a, b, 1, 'Time(sec)', th_label,
                      'Time History', 'time_history', 'PSD')

    plt.gca().set_autoscale_on(False)

    directory = "PSD/magnitude"

    if not os.path.exists(directory):
        os.makedirs(directory)

    ave_line = [np.mean(full)]*len(freq)
    thresh = (np.mean(full))*1.5

    peaks, _ = find_peaks(full, height=thresh)
    plt.figure(2, figsize=[12, 8])
    plt.plot(freq, full)
    plt.plot(freq[peaks], full[peaks], "x")
    plt.plot(freq, ave_line)

    plt.plot(freq[idx], full[idx], '*')
    title_string = 'PSD ' + \
        str("%6.3g" % rms)+' GRMS Overall ' + str(datetime.datetime.now())
    plt.title(title_string)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(psd_label)
    plt.xlabel(' Frequency (Hz) ')
    plt.grid(True)
    timestamp = str(datetime.datetime.now())
    timestamp = timestamp.replace(' ', '-')
    timestamp = timestamp.replace('.', '_')
    timestamp = timestamp.replace(':', '-')
    stitle = 'power_spectral_density' + timestamp + ".png"
    # savepath = directory + "/" + stitle
    # stitle = 'power_spectral_density' + str(datetime.datetime.now()) + ".png"
    savepath = os.path.join(directory, stitle)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(savepath)
    plt.show(block=False)
    plt.pause(1)
    plt.close('all')

########################################################################


def psd_post(freq, full, rms, run_psd, output_psd):
    print(" ")
    print(" Overall RMS = %10.3g " % rms)
    print(" Three Sigma = %10.3g " % (3*rms))

    idx = argmax(full)

    print(" ")
    print(" Maximum:  Freq=%8.4g Hz   Amp=%8.4g " % (freq[idx], full[idx]))

    # print(" ")
    # print(" Write PSD data to file? 1=yes 2=no")
    # iacc = GetInteger2()
    iacc = 1

    if(iacc == 1):

        if(run_psd == 1):

            # print(" ")
            # print(" Find output dialog box")

            # root = tk.Tk()
            # root.withdraw()
            # output_file_path = asksaveasfilename(
            #     parent=root, title="Enter the PSD output filename: ")
            # output_file = output_file_path.rstrip('\n')
            # output_file = "psd_data"
            output_file = input(
                "Enter file name(without extension) to store psd data : ")
            output_file = 'PSD/' + output_file + '.csv'
            mH = len(freq)
            WriteData2(mH, freq, full, output_file)
            run_psd = 0
        else:
            output_file = output_psd
            mH = len(freq)
            WriteData2(mH, freq, full, output_file)

    return (run_psd, output_file, idx)
########################################################################
