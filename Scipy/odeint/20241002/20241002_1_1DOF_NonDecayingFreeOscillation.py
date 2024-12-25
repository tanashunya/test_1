'''
1自由度非減衰自由振動
'''

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint


def f(var, t, m, k):
    ''' Harmonic oscillator '''

    # x, vに分離
    x = var[0]
    v = var[1]

    # 運動方程式
    a = -(k / m) * x

    # 結果
    output = np.array([v, a])

    return output

def plot(t, x):
    ''' Plot '''

    # フォントの種類とサイズを設定する。
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする。
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右に目盛線を付ける。
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')

    # 軸のラベルを設定する。
    ax1.set_xlabel('t[s]')
    ax1.set_ylabel('x[m]')

    # データプロットの準備とともに、ラベルと線の太さ、凡例の設置を行う。
    ax1.plot(t, x, label='Displacement', lw=1, color='red')

    # レイアウト設定
    fig.tight_layout()

    # グラフを表示する。
    plt.show()
    plt.close()

if __name__ == '__main__':
    ''' Main '''

    # モデルを定義
    m = 1.0
    k = 10.0

    # 解析時間
    dt = 1e-3
    t_max = 10.0
    t = np.arange(0.0, t_max, dt)

    # 初期条件
    x0 = 1.0
    v0 = 0.0
    var = [x0, v0]

    # 微分方程式の近似解法（過渡応答）
    var = odeint(f, var, t, args=(m, k))

    # 結果を分離
    x = var.T[0]
    v = var.T[1]

    # プロット
    plot(t, x)