'''
1次元減衰自由振動
'''

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

def f(var, t, m, k, c):
    '''Harmonic oscillator'''
    x = var[0]
    v = var[1]

    # 運動方程式
    a = -(k / m) * x - (c /m) * v

    # 結果(1D配列にする必要があり)
    output = np.array([v, a])

    return output

def plot(t, x):
    '''plot'''

    # フォントの種類とサイズを設定する
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右にメモリ線をつける
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # 軸のラベルを指定する
    ax1.set_xlabel('t[s]')
    ax1.set_ylabel('x[m]')

    ax1.plot(t, x, label='Displacement', lw=1, color='red')

    fig.tight_layout()

    plt.show()
    plt.close()

if __name__ == '__main__':
    # 定数
    m = 1.0
    k = 10.0
    c = 1.0

    # 解析時間
    dt = 1e-3
    t_max = 10.0
    t = np.arange(0.0, t_max, dt)

    # 初期条件
    x0 = 1.0
    v0 = 0.0
    var = [x0, v0]

    # odeint
    var = odeint(f, var, t, args=(m, k, c))

    # 結果を分離
    x = var.T[0]
    v = var.T[1]

    # プロット
    plot(t, x)


