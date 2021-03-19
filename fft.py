import serial
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft
import tkinter as tk
from tkinter import filedialog
import time
import csv
import pandas as pd
import sys
from time import sleep
import urllib.request
import blynklib
from scipy.fft import fft, fftfreq
from scipy.signal import hann, blackman, hamming, square
import traceback
import matplotlib.animation as animation
from matplotlib import style
import psd_lib
import pathlib
import tompy
# from sys import stdinpip


def select_window():
    print("For Hann - 1.")
    print("For Rectangular - 2.")
    print("For Blackman - 3.")
    print("For Hamming - 4")
    print("None - any other key.")
    select = int(input("Enter"))
    if select == 1:
        window = "Hann"
    elif select == 2:
        window = "Rectangular"
    elif select == 3:
        window = "Blackman"
    elif select == 4:
        window = "Hamming"
    else:
        window = None
    return (select, window)


def get_user_input(message):
    value = int(input(message))
    if value == 1:
        return True
    return False


def get_data(url):
    """Returns data"""
    n = urllib.request.urlopen(url).read()
    # get the raw html data in bytes (sends request and warn our esp8266)
    data = n.decode("utf-8")  # convert raw html bytes format to string :3
    return data


def max_yf_check(num_samples, yf_ax, T, N):
    """Creates a max_frequency.csv if fft is calculated for the first time, otherwise finds the positional maximum value and replaces the csv."""
    dirname = pathlib.Path(__file__).parent
    filename = dirname / 'max_frequency.csv'
    
    if not filename.is_file():
        max_yf_ax = np.zeros(num_samples//2)+0j
        # max_yf_ay = np.zeros(num_samples)+0j
        # max_yf_az = np.zeros(num_samples)+0j
    
    else:
        with open(filename) as filename_csv:
            df = pd.DataFrame(list(csv.reader(filename_csv)), columns = ['Amplitude', 'Frequency'])
            max_yf_ax = df['Amplitude'].values.astype(complex)

    # Comparing current yf with max_yf
    

    for i in range(len(max_yf_ax)):
        if abs(yf_ax[i]) > abs(max_yf_ax[i]):
            max_yf_ax[i] = yf_ax[i]

    # for i in range(len(max_yf_ay)):
        # if abs(yf_ay[i]) > abs(max_yf_ay[i]):
            # max_yf_ay[i] = yf_ay[i]

    # for i in range(len(max_yf_az)):
        # if abs(yf_az[i]) > abs(max_yf_az[i]):
            # max_yf_az[i] = yf_az[i]

    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    df = pd.DataFrame(columns = ['Amplitude', 'Frequency'])
    df['Amplitude'] = max_yf_ax
    df['Frequency'] = xf 
    
    
    
    # if filename.is_file() and (not (yf_ax == max_yf_ax).all()):
    # print('Creating or updating max_frequency.csv')
    df.to_csv('max_frequency.csv', header=False, index=False)
    
        
def save_to_csv(df, filename):
    """Saves a dataframe to a file with the given filename. If file already exists, the dataframe will be appended."""
    dirname = pathlib.Path(__file__).parent
    file = dirname / filename
    
    if not file.is_file():
        # print('Creating {}'.format(filename))
        df.to_csv(filename, header=False, index=False)
    else:
        # print('Appending {}'.format(filename))
        df.to_csv(filename, mode='a', header=False, index=False) 
        

def calculate_save_fft(df_ax, delta_time, num_samples, window):
    """Calculates fft of a dataframe, saves the calculations to a csv. Creates a max_frequency.csv if fft is calculated for the first time, otherwise finds the positional maximum value and replaces the csv."""
    root = tk.Tk()
    root.withdraw()

    t = df_ax['Time']
    # delta_time = (t[3]-t[2]+t[2]-t[1] + t[4]-t[3] + t[5]-t[4]+ t[6]-t[5])/5000
    # print(delta_time)

    x_ax = df_ax['Accelaration']
    # x_ay = df_ay['Accelaration']
    # x_az = df_az['Accelaration']

    # Determine variables
    N = int(np.prod(t.shape))  # length of the array
    Fs = 1/delta_time  # sample rate (Hz)
    T = 1/Fs


    # Compute FFT
    if window:
        if window == "Hann":
            w = hann(num_samples)
        if window == "Rectangular":
            w = square(num_samples)
        if window == "Blackman":
            w = blackman(num_samples)
        if window == "Hamming":
            w = hamming(num_samples)
        yf_ax = fft(x_ax.values*w)

    else:
        yf_ax = fft(x_ax)


    yf_ax = yf_ax [1:(N//2) + 1]
    xf = np.linspace(0.0, 1.0/(2.0*delta_time), N//2)
    # yf_ay = fft(x_ay)
    # yf_az = fft(x_az)

    # yf_ax[0] = 0
    # yf_ay[0] = 0
    # yf_az[0] = 0

    # yf_ax[1] = 0
    # yf_ay[1] = 0
    # yf_az[1] = 0

    # yf_ax[-1] = 0
    # yf_ay[-1] = 0
    # yf_az[-1] = 0

    # df_ax['Freq'] = xf
    # df_ay['Freq'] = xf
    # df_az['Freq'] = xf
    
    df_freq = pd.DataFrame(columns=['Amplitude', 'Frequency'])
    df_freq['Amplitude'] = yf_ax
    df_freq['Frequency'] = xf
    
    if window:
        save_to_csv(df_freq, 'fft {}.csv'.format(window))
    else:
        save_to_csv(df_freq, 'fft.csv')
    max_yf_check(num_samples, yf_ax, T, N)
    

def psd_select_options():
    
    # print(" ")
    # print(" Remove mean:  1=yes  2=no ")

    # mr_choice = tompy.GetInteger2()
    mr_choice = 1

    # print(" ")
    # print(" Select Window: 1=Rectangular 2=Hanning ")

    # h_choice = tompy.GetInteger2()
    h_choice = 1
    
    # print(' ')
    # print(' Choose the Number of Segments:  ')
    # s = stdin.readline()
    s = 1

    return mr_choice, h_choice, s


def create_dataframes(delta_time, input_time, url):
    """Returns a dataframe with acceleration in 'Accelaration' column and time in 'Time' column,"""

    for i in range(2):
        get_data(url)
    
    first = True
    arduinoData = get_data(url)
    arduinoData = arduinoData.split()
    df_ax = pd.DataFrame(columns=['Time', 'Accelaration'])
    df_ax['Accelaration'] = list(map(int, arduinoData[::2])) 
    df_ax['Time'] = list(map(float, arduinoData[1::2])) 
    df_ax['Time'] = df_ax['Time'] / 1000000.0
    num_samples = len(df_ax['Time'])

    # if arduinoData == b'' or arduinoData == b'0\n':
    #     b_blank += 1
    #     if b_blank == 10:
    #         sys.exit('Blank bytes.')
    # else:
    #     b_blank = 0

    # if(arduinoData == b'MPU6050 connection failed\r\n'):
    #     sys.exit(arduinoData)

    # <optional> split datas we got. (if you programmed it to send more than one value) It splits them into seperate list elements.
    
    # arduinoData = arduinoData.split()
    # arduinoData_ax = int(arduinoData[0])
    # arduinoData_time = int(arduinoData[1])
    # arduinoData_az = int(arduinoData[2])
    
    # if first or arduinoData == b'':
    #     df_ax = df_ax.append({'Time': arduinoData_time, 'Accelaration': 0}, ignore_index=True)
        # df_ay = df_ay.append({'Time':input_time, 'Accelaration':0}, ignore_index=True)
        # df_az = df_az.append({'Time':input_time, 'Accelaration':0}, ignore_index=True)
    #     if first:
    #         first = False
    # else:
        
    #     df_ax = df_ax.append(
    #         {'Time': arduinoData_time, 'Accelaration': int(arduinoData_ax)}, ignore_index=True)
        # df_ay = df_ay.append({'Time':input_time, 'Accelaration':int(arduinoData_ay)}, ignore_index=True)
        # df_az = df_az.append({'Time':input_time, 'Accelaration':int(arduinoData_az)}, ignore_index=True)

    return (df_ax, num_samples)

    
    #     input_time = round(input_time + delta_time, 3)
        # x = df_ax['Time']
        # y = df_ax['Accelaration']
    

def calculate_save_psd(a,b,num, mr_choice, h_choice, s):
    print("The input file must have two columns: time(sec) & amplitude.")

    # Calculate PSD
    # a = Time in list format
    # b = Acceleration in list format
    # num = Number of data points for a single calculation

    sr,dt,dur=psd_lib.READ_DATA().read_and_stats_custom(a,b,num)
    mmm,NW,df,mH,mr_choice,h_choice=psd_lib.select_options(num, dt, mr_choice, h_choice, s)  
    rms,freq,full=psd_lib.PSD(mmm,NW,mH,df,b,mr_choice,h_choice).psd_core()
    # Save PSD Data
    psd_lib.psd_post(freq,full,rms)
    return freq,full,rms, a, b

def plot_psd(freq,full,rms, a, b):
    '''Plot PSD'''
    # GET N VALUES FROM psd_data.csv
    psd_lib.psd_plots(freq,full,rms, a, b)
    
# def calculate_plot_rms(*args):
    # tic = time.process_time()
    # w = np.int(np.floor(Fs)); #width of the window for computing RMS
    # steps = np.int_(np.floor(N/w)); #Number of steps for RMS
    # t_RMS = np.zeros((steps,1)); #Create array for RMS time values
    # x_RMS = np.zeros((steps,1)); #Create array for RMS values
    # for i in range (0, steps):
    #     t_RMS[i] = np.mean(t[(i*w):((i+1)*w)]);
    #     x_RMS[i] = np.sqrt(np.mean(x[(i*w):((i+1)*w)]**2));  
    # plt.figure(2)  
    # plt.plot(t_RMS, x_RMS)
    # plt.xlabel('Time (seconds)')
    # plt.ylabel('RMS Accel (g)')
    # plt.title('RMS - ' + file_path)
    # plt.grid()
    # toc = time.process_time()
    # print("RMS Time:",toc-tic)




def plot_fft(window, num_samples):
    """Plot FFT."""
    if window:
        df = pd.read_csv("fft {}.csv".format(window))
    else:
        df = pd.read_csv("fft.csv")
    df = df.tail(num_samples)
    df.iloc[:,0] = df.iloc[:,0].astype(complex)
    df.iloc[:,1] = df.iloc[:,1].astype(float)
    yf_ax = df.iloc[:,0] # Amplitude
    xf = df.iloc[:,1] # Frequency
    N = num_samples//2 # Already depricated the mirrored half values before saving into CSV.
    plt.figure(num=None, figsize=(20, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(xf[:int(N)], 2.0/N * np.abs(yf_ax[0:int(N)]), '-g')
    
    if window:
        plt.legend(['FFT windowed with ' + window, 'FFT'])
    else:
        plt.legend(['FFT'])

    plt.grid()
    plt.xlabel('Frequency of ax(Hz)')
    plt.ylabel('Accel (g)')
    if window:
        plt.title('FFT with {} window'.format(window))
    else:
        plt.title('FFT')
    plt.show()
    

def choices():
    choice = int(input("Plots you want to check - \n1. FFT\n2. PSD\n"))
    if choice == 1:
        return 'fft'
    return 'psd'

if __name__ == '__main__':    
    # ESP's IP, ex: http://192.168.102/ (Check serial console while uploading the ESP code, the IP will be printed)
    url = "http://192.168.143.189"
    delta_time = .0017
    input_time = 0
    select, window = select_window()
    mr_choice, h_choice, s = psd_select_options()
    choice = choices()
    num_samples = 1024
    while True:
        df_ax, num_samples = create_dataframes(delta_time, input_time, url) 
        save_to_csv(df_ax, 'acceleration_data.csv') # Saving acceleration data
        calculate_save_fft(df_ax, delta_time, num_samples, window)
        if choice == 'fft':
            plot_fft(window, num_samples)
        freq,full,rms, a, b = calculate_save_psd(df_ax['Time'], df_ax['Accelaration'], num_samples, mr_choice, h_choice, s)
        if choice == 'psd':
            plot_psd(freq,full,rms, a, b)



# TODO
# 4. Third script - Plots PSD
# 5. RMS integrate
# 6. Check if required - Delta calculate from initial time values and modify the time data accordingly.
# Change file once a particular limit is reached
# Skip data leakage - twice get_data(url) wasted
# PSD Select options function in fft.py correct - Manual input
#  Is the input data dimension Accel(G) ? - Correct - Manual input