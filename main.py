from p_psd import *
from p_fft import *
import urllib.request
import pathlib
import pandas as pd


def choices():
    choice = int(input("Plots you want to check - \nFFT : 1\nPSD : 2\n"))
    if choice == 1:
        return 'fft'
    else:
        return 'psd'


def psd_main():
    print(" The input file must have two columns: time(sec) & amplitude ")

    a, b, num, sr, dt, dur = READ_DATA().read_and_stats()

    mmm, NW, df, mH, mr_choice, h_choice = select_options(num, dt)

    rms, freq, full = PSD(mmm, NW, mH, df, b, mr_choice,
                          h_choice).psd_core()

    psd_post(freq, full, rms)

    psd_plots(a, b, freq, full, rms)


def main_fft():
    print(" ")

    print("The file must have two columns: time(sec) & amplitude")

    a, b, num, sr, dt, dur = READ_DATA().read_and_stats()

    print(" ")
    print(" Remove mean:  1=yes  2=no ")

    imr = GetInteger2()

    ########################################################################

    idx, freq, ff, z, zz, ph, nhalf, df, num_fft = FFT(a, b, imr).fft_data()

    ########################################################################

    print(" ")
    print(" Maximum:  Freq=%8.4g Hz   Amp=%8.4g " % (ff[idx], zz[idx]))

    print(" ")
    print("  Export data files:  1=yes 2=no ")

    idf = GetInteger2()

    if(idf == 1):

        print(" ")
        print(" Find output dialog box")

        root = tk.Tk()
        root.withdraw()
        output_file_path = asksaveasfilename(
            parent=root, title="Enter the output FFT (freq, real, imag): ")
        output_file = output_file_path.rstrip('\n')
        WriteData3(num_fft, freq, z.real, z.imag, output_file)

        root = tk.Tk()
        root.withdraw()
        output_file_path = asksaveasfilename(
            parent=root, title="Enter the output FFT filename (freq, mag, phase(rad)): ")
        output_file = output_file_path.rstrip('\n')
        WriteData3(nhalf, ff, zz, ph, output_file)

    print(" ")
    print(" view plots ")

    time_history_plot(a, b, 1, 'Time(sec)', 'Amplitude',
                      'Time History', 'time_history')

    plt.figure(2)
    plt.plot(ff, zz)
    plt.grid(True)
    plt.title(' FFT Magnitude ')
    plt.ylabel(' Amplitude ')
    plt.xlabel(' Frequency (Hz) ')
    plt.grid(True, which="both")
    plt.savefig('fourier_magnitude')
    plt.draw()

    plt.figure(3)
    plt.plot(ff, ph*(180./pi))
    plt.grid(True)
    plt.title(' FFT Phase ')
    plt.ylabel(' Phase (deg) ')
    plt.xlabel(' Frequency (Hz) ')
    plt.grid(True, which="both")
    plt.savefig('fourier_phase')
    plt.draw()

    plt.show()


# def select_window():
#     print("For Hann - 1.")
#     print("For Rectangular - 2.")
#     print("For Blackman - 3.")
#     print("For Hamming - 4")
#     print("None - any other key.")
#     select = int(input("Enter"))
#     if select == 1:
#         window = "Hann"
#     elif select == 2:
#         window = "Rectangular"
#     elif select == 3:
#         window = "Blackman"
#     elif select == 4:
#         window = "Hamming"
#     else:
#         window = None
#     return (select, window)


def get_data(url):
    n = urllib.request.urlopen(url).read()
    data = n.decode("utf-8")
    return data


def create_dataframes(url):

    arduinoData = get_data(url)
    arduinoData = arduinoData.split()
    df_ax = pd.DataFrame(columns=['Time', 'Accelaration'])
    df_ax['Accelaration'] = list(map(int, arduinoData[::2]))
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


url = "http://192.168.143.189"
choice = choices()
while True:
    df_ax, num_samples = create_dataframes(url)
    save_to_csv(df_ax, 'acceleration_data.csv')  # Saving acceleration data
    if(choice == 'fft'):
        main_fft()
    else:
        psd_main()
