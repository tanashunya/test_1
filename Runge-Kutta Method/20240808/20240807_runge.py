import numpy as np
from scipy import fftpack
from matplotlib import pyplot as plt
from scipy.integrate import odeint

def model():
    ''' Define model '''

    # 質量・減衰・剛性の集中定数を設定する
    m1 = 0.5
    m2 = 0.3
    c1 = 0
    c2 = 0
    c12 = 0
    k1 = 10000
    k2 = 10000
    k12 = 45

    M = np.array([[m1, 0],
                  [0, m2]])
    C = np.array([[c1 + c12, -c12    ],
                  [-c12    , c2 + c12]])
    K = np.array([[k1 + k12, -k12    ],
                  [-k12    , k2 + k12]])

    return M, K, C

def f(var, t, M, K, C, A, f0):
    ''' Harmonic oscillator '''

    # 質量行列の逆行列
    M_inv = np.linalg.inv(M)

    # 外力
    F = np.zeros(len(M))
    F[-1] = A * np.sin(2 * np.pi * f0)

    # x, vに分離
    x = var[0:int(len(var) / 2)]
    v = var[int(len(var) / 2)::]

    # 運動方程式
    a = - (M_inv @ K) @ x - (M_inv @ C) @ v + (M_inv @ F)

    # 結果を平坦化（odeint用に1Dにする）
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

    return amp, freq

def plot(t, x, force, amp, freq):
    ''' プロット '''

    # フォントの種類とサイズを設定する。
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする。
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右に目盛線を付ける。
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(221)
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax2 = fig.add_subplot(223)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax3 = fig.add_subplot(122)
    ax3.yaxis.set_ticks_position('both')
    ax3.xaxis.set_ticks_position('both')

    # 軸のラベルを設定する。
    ax1.set_xlabel('t[s]')
    ax1.set_ylabel('x[m]')
    ax2.set_xlabel('t[s]')
    ax2.set_ylabel('F[N]')
    ax3.set_xlabel('Driving frequency[Hz]')
    ax3.set_ylabel('Compliance[m/N]')

    # スケールの設定をする。
    ax3.set_xlim(0.0, 10.0)
    ax3.set_yscale('log')

    # データプロットの準備とともに、ラベルと線の太さ、凡例の設置を行う。
    ax1.plot(t, x, label='Displacement', lw=1, color='red')
    ax2.plot(t, force, label='Force', lw=1, color='red')
    ax3.plot(freq, amp, label='FRF', lw=1, color='red')

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
    Amplitude = 1.0
    f0 = 1.0
    force = Amplitude * np.sin(2 * np.pi * f0)

    # 微分方程式の近似解法（過渡応答）
    var = odeint(f, var, t, args=(M, K, C, Amplitude, f0))

    # 結果を分離（odeint用の1D配列から物理量毎に分離）
    x = var.T[0:int(len(var.T) / 2)]
    v = var.T[int(len(var.T) / 2)::]

    # amp, freq_axis = frf(force, x[-1], 1 / dt)

    # プロット
    # plot(t, x[-1], force, amp, freq_axis)