import datetime
import matplotlib.pyplot as plt
from scipy.signal import windows
import pandas as pd
import pathlib
import urllib.request
from p_fft import *
from p_psd import *
import os


def choices():
    choice = int(
        input("Plots you want to check - \nFFT : 1\nPSD : 2\nEnter : "))
    if choice == 1:
        return 'fft'
    else:
        return 'psd'


def psd_main(window, run_psd, output_psd):
    # print(" The input file must have two columns: time(sec) & amplitude ")
    directory = "PSD"

    if not os.path.exists(directory):
        os.makedirs(directory)

    a, b, num, sr, dt, dur = READ_DATA_PSD().read_and_stats()

    mmm, NW, df, mH, mr_choice, h_choice = select_options(num, dt, window)

    rms, freq, full = PSD(mmm, NW, mH, df, b, mr_choice,
                          h_choice).psd_core()

    run_psd, output_psd, idx = psd_post(freq, full, rms, run_psd, output_psd)

    psd_plots(a, b, freq, full, rms, idx)

    return (run_psd, output_psd)


def fft_main(window, run_fft, output_file_file1, output_file_file2):

    directory = "FFT"

    if not os.path.exists(directory):
        os.makedirs(directory)

    # print(" ")

    # print("The file must have two columns: time(sec) & amplitude")

    a, b, num, sr, dt, dur, ave, rms = READ_DATA_FFT().read_and_stats()
    # print(type(a))

    # print(" ")
    # print(" Remove mean:  1=yes  2=no ")

    # imr = GetInteger2()
    imr = 1

    if window == "Hann":
        w = windows.hann(num)
    elif window == "Rectangular":
        w = windows.boxcar(num)
    elif window == "Blackman":
        w = windows.blackman(num)
    elif window == "Hamming":
        w = windows.hamming(num)
    else:
        w = 1

    b = [b[i] * w[i] for i in range(len(b))]
    # a = [x * w for x in a]

    ########################################################################

    idx, freq, ff, z, zz, ph, nhalf, df, num_fft = FFT(a, b, imr).fft_data()

    ########################################################################
    # max_y = max(zz)  # Find the maximum y value
    # # Find the x value corresponding to the maximum y value
    # max_x = ff[zz.argmax()]
    # print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    # print(max_x)
    # print(max_y)

    data_to = str(" Maximum:  Freq : %8.4g Hz   Amp : %8.4g\nRMS : %8.4g   Mean : %8.4g" % (
        ff[idx], zz[idx], rms, ave))
    print(" ")
    # print(" Maximum:  Freq=%8.4g Hz   Amp=%8.4g " % (ff[idx], zz[idx]))
    print(data_to)

    # print(" ")
    # print("  Export data files:  1=yes 2=no ")

    # idf = GetInteger2()
    idf = 1

    if(idf == 1):

        # print(" ")
        # print(" Find output dialog box")

        if(run_fft == 1):
            # root = tk.Tk()
            # root.withdraw()
            # output_file_path = asksaveasfilename(
            #     parent=root, title="Enter the output FFT (freq, real, imag): ")
            # output_file = output_file_path.rstrip('\n')
            # output_file = "fft_data_file1"
            output_file_file1 = input(
                "Enter file name(without extension) to store frequency and complex data : ")
            output_file_file1 = 'FFT/' + output_file_file1 + '.csv'
            WriteData3(num_fft, freq, z.real, z.imag, output_file_file1)
            run_fft = 0

            # root = tk.Tk()
            # root.withdraw()
            # output_file_path = asksaveasfilename(
            #     parent=root, title="Enter the output FFT filename (freq, mag, phase(rad)): ")
            # output_file = output_file_path.rstrip('\n')
            # output_file = "fft_data_file2"
            output_file_file2 = input(
                "Enter file name(without extension) to store frequency and amplitude data : ")
            output_file_file2 = 'FFT/' + output_file_file2 + '.csv'
            WriteData3(nhalf, ff, zz, ph, output_file_file2)
        else:
            WriteData3(num_fft, freq, z.real, z.imag, output_file_file1)
            WriteData3(nhalf, ff, zz, ph, output_file_file2)

    print(" ")
    print(" view plots ")

    time_history_plot(a, b, 1, 'Time(sec)', 'Amplitude',
                      'Time History', 'time_history', 'FFT')

    directory = "FFT/magnitude"

    if not os.path.exists(directory):
        os.makedirs(directory)

    plt.figure(2, figsize=[12, 8])
    plt.plot(ff, zz)
    plt.plot(ff[idx], zz[idx], '*')
    plt.grid(True)
    ptitle = 'FFT Magnitude ' + str(datetime.datetime.now())
    plt.title(ptitle)
    plt.ylabel(' Amplitude ')
    plt.xlabel(' Frequency (Hz) \n' + data_to)
    plt.grid(True, which="both")
    timestamp = str(datetime.datetime.now())
    timestamp = timestamp.replace(' ', '-')
    timestamp = timestamp.replace('.', '_')
    timestamp = timestamp.replace(':', '-')
    stitle = 'fourier_magnitude' + timestamp + ".png"
    savepath = os.path.join(directory, stitle)
    plt.savefig(savepath)
    plt.draw()

    directory = "FFT/phase"

    if not os.path.exists(directory):
        os.makedirs(directory)

    plt.figure(3, figsize=[12, 8])
    plt.plot(ff, ph*(180./pi))
    plt.grid(True)
    ptitle = 'FFT Phase ' + str(datetime.datetime.now())
    plt.title(ptitle)
    plt.ylabel(' Phase (deg) ')
    plt.xlabel(' Frequency (Hz) ')
    plt.grid(True, which="both")
    timestamp = str(datetime.datetime.now())
    timestamp = timestamp.replace(' ', '-')
    timestamp = timestamp.replace('.', '_')
    timestamp = timestamp.replace(':', '-')
    stitle = 'fourier_phase' + timestamp + ".png"
    savepath = os.path.join(directory, stitle)
    plt.savefig(savepath)
    plt.draw()

    plt.show(block=False)
    plt.pause(1)
    plt.close('all')
    return (run_fft, output_file_file1, output_file_file2)


def select_window():
    print("For Rectangular : 1")
    print("For Hann : 2")
    print("For Blackman : 3")
    print("For Hamming : 4")
    print("None - any other key.")
    select = int(input("Enter : "))
    if select == 1:
        window = "Rectangular"
    elif select == 2:
        window = "Hann"
    elif select == 3:
        window = "Blackman"
    elif select == 4:
        window = "Hamming"
    else:
        window = None
    return (window)


def get_data(url):
    n = urllib.request.urlopen(url).read()
    data = n.decode("utf-8")
    return data


# def create_dataframes(url):
#     infile = open("arduino.txt", "r")
#     arduinoData = infile.readline()
#     arduinoData = [float(idx) for idx in arduinoData.split(' ')]
#     df_ax = pd.DataFrame(columns=['Time', 'Accelaration'])
#     df_ax['Accelaration'] = list(map(float, arduinoData[:: 2]))
#     df_ax['Time'] = list(map(float, arduinoData[1:: 2]))
#     df_ax['Time'] = df_ax['Time'] / 1000000.0
#     num_samples = len(df_ax['Time'])

#     return (df_ax, num_samples)


def create_dataframes(url):

    arduinoData = get_data(url)
    arduinoData = arduinoData.split()
    df_ax = pd.DataFrame(columns=['Time', 'Accelaration'])
    df_ax['Accelaration'] = list(map(float, arduinoData[::2]))
    df_ax['Time'] = list(map(float, arduinoData[1::2]))
    df_ax['Time'] = df_ax['Time'] / 1000000.0
    num_samples = len(df_ax['Time'])

    return (df_ax, num_samples)


def save_to_csv(df, filename):
    dirname = pathlib.Path(__file__).parent
    file = dirname / filename
    if not file.is_file():
        df.to_csv(filename, header=False, index=False)
    else:
        df.to_csv(filename, mode='a', header=False, index=False)


run_fft = 1
run_psd = 1
output_file_file1 = None
output_file_file2 = None
output_psd = None

url = "http://192.168.143.189"
window = select_window()
choice = choices()
while True:
    df_ax, num_samples = create_dataframes(url)
    save_to_csv(df_ax, 'acceleration_data.csv')  # Saving acceleration data
    if(choice == 'fft'):
        run_fft, output_file_file1, output_file_file2 = fft_main(
            window, run_fft, output_file_file1, output_file_file2)
    else:
        run_psd, output_psd = psd_main(
            window, run_psd, output_psd)
