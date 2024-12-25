
import numpy as np
from matplotlib import pyplot as plt
import scipy
import scipy.fftpack
from scipy.integrate import odeint

def f(var, t, m, k, c, A, alpha, f0):
    '''Harmonic oscillator'''
    x = var[0]
    v = var[1]

    # 外力を計算
    force = A * np.sin(2 * np.pi * (f0 * t + (alpha / 2) * t ** 2))

    # 運動方程式
    a = -(k / m) * x - (c /m) * v + force

    # 結果(1D配列にする必要があり)
    output = np.array([v, a])

    return output

def frf(input, output, samplerate):
    '''周波数応答関数（FRF）'''
    # 入出力信号のフーリエ変換
    fft_input = scipy.fftpack.fft(input)
    fft_output = scipy.fftpack.fft(output)

    # FRFを計算
    h_io = (fft_output * fft_input.conjugate()) / (fft_input * fft_input.conjugate())
    amp = np.sqrt((h_io.real ** 2) + (h_io.imag ** 2)) / (len(input) / 2)
    phase = np.arctan2(h_io.imag, h_io.real)
    phase = np.degrees(phase)
    freq = np.linspace(0, samplerate, len(input))

    return amp, freq

def plot(t, x, force, amp, freq):
    '''plot'''

    # フォントの種類とサイズを設定する
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右にメモリ線をつける
    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(221)
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    
    ax2 = fig.add_subplot(223)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')

    ax3 = fig.add_subplot(122)
    ax3.yaxis.set_ticks_position('both')
    ax3.xaxis.set_ticks_position('both')

    # 軸のラベルを指定する
    ax1.set_xlabel('t[s]')
    ax1.set_ylabel('x[m]')

    ax2.set_xlabel('t[s]')
    ax2.set_ylabel('F[N]')

    ax3.set_xlabel('Driving frequency[Hz]')
    ax3.set_ylabel('Compliance[m/N]')

    # スケールを設定する
    ax3.set_xlim(0.0, 50.0)
    ax3.set_yscale('log')

    ax1.plot(t, x, label='Displacement', lw=1, color='red')
    ax2.plot(t, force, label='Force', lw=1, color='red')
    ax3.plot(freq, amp, label='FRF', lw=1, color='red')
    
    # レイアウト設定
    fig.tight_layout()

    plt.show()
    plt.close()

if __name__ == '__main__':
    # 定数
    m = 1.0
    k = 1000.0
    c = 1.0

    # 解析時間
    dt = 1e-3
    t_max = 20.0
    t = np.arange(0.0, t_max, dt)

    # 初期条件
    x0 = 0.0
    v0 = 0.0
    var = [x0, v0]

    # 外力の計算
    f0 = 1
    f1 = 50
    A = 100.0
    alpha = (f1 - f0) / t_max
    force = A * np.sin(2 * np.pi * (f0 * t + (alpha / 2) * t ** 2))

    # odeint
    var = odeint(f, var, t, args=(m, k, c, A, alpha, f0))

    # 結果を分離
    x = var.T[0]
    v = var.T[1]

    # FRFを計算
    amp, freq = frf(force, x, 1/dt)

    # プロット
    plot(t, x, force, amp, freq)


