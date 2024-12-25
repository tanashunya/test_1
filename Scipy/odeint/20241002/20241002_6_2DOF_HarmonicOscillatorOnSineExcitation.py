import numpy as np
from scipy import fftpack
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

def model():
    ''' Define model '''

    # 質量・減衰・剛性の集中定数を設定する
    m1 = 100
    m2 = 50
    k1 = 5.0e4
    k2 = 1.0e4
    c1 = 5
    c2 = 5

    M = np.array([[m1, 0],
                  [0, m2]])
    K = np.array([[k1 + k2, -k2],
                  [-k2, k2]])
    C = np.array([[c1 * c2, -c2],
                  [-c2, c2]])

    return M, K, C

def f(t, var, M, K, C, A, k, f0):
    ''' Harmonic oscillator '''

    # 質量行列の逆行列
    M_inv = np.linalg.inv(M)

    # 外力
    F = np.zeros(len(M))
    F[-1] = A * np.sin(2 * np.pi * (f0 * t + (k / 2) * t ** 2))

    # x, vに分離
    x = var[0:int(len(var) / 2)]
    v = var[int(len(var) / 2)::]

    # 運動方程式
    a = - (M_inv @ K) @ x - (M_inv @ C) @ v + (M_inv @ F)

    # 結果を平坦化（solve_ivp用に1Dにする）
    output = np.concatenate([v, a])

    return output

def frf(input, output, samplerate):
    ''' 周波数応答関数(FRF) '''

    # 入出力信号のフーリエ変換
    fft_i = fftpack.fft(input)
    fft_o = fftpack.fft(output)

    # FRFを計算
    h_io = (fft_o * fft_i.conjugate()) / (fft_i * fft_i.conjugate())
    amp = np.sqrt((h_io.real ** 2) + (h_io.imag ** 2))
    amp = amp / (len(input) / 2)
    phase = np.arctan2(h_io.imag, h_io.real)
    phase = np.degrees(phase)
    freq = np.linspace(0, samplerate, len(input))

    return amp, freq, fft_o

def plot(t, x, force, amp, freq, fft_i):
    ''' プロット '''

    # フォントの種類とサイズを設定する。
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする。
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右に目盛線を付ける。
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(221)
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax2 = fig.add_subplot(223)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax3 = fig.add_subplot(222)
    ax3.yaxis.set_ticks_position('both')
    ax3.xaxis.set_ticks_position('both')
    ax4 = fig.add_subplot(224)
    ax4.yaxis.set_ticks_position('both')
    ax4.xaxis.set_ticks_position('both')

    # 軸のラベルを設定する。
    ax1.set_xlabel('t[s]')
    ax1.set_ylabel('x[m]')
    ax2.set_xlabel('t[s]')
    ax2.set_ylabel('F[N]')
    ax3.set_xlabel('Frequency[Hz]')
    ax3.set_ylabel('FFT Amplitude')
    ax4.set_xlabel('Driving frequency[Hz]')
    ax4.set_ylabel('Compliance[m/N]')

    # スケールの設定をする。
    ax3.set_xlim(0.0, 15)

    ax4.set_xlim(0.0, 10.0)
    ax4.set_yscale('log')

    # データプロットの準備とともに、ラベルと線の太さ、凡例の設置を行う。
    ax1.plot(t, x, label='Displacement', lw=0.05, color='red')
    ax2.plot(t, force, label='Force', lw=0.05, color='red')
    ax3.plot(freq, np.abs(fft_i), label='FFT', lw=1, color='red')
    ax4.plot(freq, amp, label='FRF', lw=1, color='red')

    # レイアウト設定
    fig.tight_layout()

    # グラフを表示する。
    plt.show()
    plt.close()

if __name__ == '__main__':
    ''' メイン処理 '''

    # モデルを定義
    M, K, C = model()

    # 解析時間
    dt = 1e-3
    t_max = 100.0
    t = np.arange(0.0, t_max, dt)

    # 初期条件
    x0 = np.zeros(len(M))
    v0 = np.zeros(len(M))
    var = np.concatenate([x0, v0])

    # 外力（forceはFRF計算用）
    A = 1.0
    f0 = 1.0
    f1 = 20.0
    k = (f1 - f0) / t_max
    force = A * np.sin(2 * np.pi * (f0 * t + (k / 2) * t ** 2))

    # 微分方程式の近似解法（過渡応答）
    sol = solve_ivp(f, [0, t_max], var, args=(M, K, C, A, k, f0), t_eval=t, method='BDF')

    # 結果を分離（solve_ivp用の1D配列から物理量毎に分離）
    x = sol.y[0:int(len(sol.y) / 2)]
    v = sol.y[int(len(sol.y) / 2)::]

    amp, freq_axis, fft_o = frf(force, x[-1], 1 / dt)

    # プロット
    plot(t, x[-1], force, amp, freq_axis, fft_o)