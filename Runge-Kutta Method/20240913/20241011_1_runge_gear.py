import numpy as np
from scipy import fftpack
from scipy.signal import chirp
import matplotlib.pyplot as plt
from RBFN_model import Fm_Linear

# 定数
m1 = 0.1    
m2 = 1.0
c1 = 0.5
c2 = 0.77
c12 = 0.0
# k1 = 200.0
# k2 = 10.0
# k_12 = 50.0
# k1 = 0.0
# k2 = 0.0
# k_12 = 0.0
k1 = 0.0
k2 = 0.0
k_12 = 0.0

M = np.array([[m1, 0],
                [0, m2]])
C = np.array([[c1 + c12, -c12],
                [-c12, c2 + c12]])
K = np.array([[k1 + k_12, -k_12],
                [-k_12, k2 + k_12]])

i_f = 0.0
i_s = 0.0

# 解析条件
t0 = 0.0
t1 = 10.0
dt = 1e-3

x0_0 = 0.0001

f0 = 4.0
f1 = 4.0
Amp = 0.1
freq = 4.0

fm_model = Fm_Linear('sum2.csv')

def Fm(z1, z2, i_f, i_s):
    # return 0.0
    fm_val = fm_model.interpolate(z1-z2, i_f, i_s)
    # print(z1-z2)
    return fm_val

def f(x, v):
    ''' 運動方程式 '''

    # 質量マトリクスの逆行列
    M_inv = np.linalg.inv(M)

    # 外力ベクトル
    F = np.zeros(len(M))
    F[0] = -m1 * x0_dotdot(t) + Fm(x[0], x[1], 0.0, 0.0)
    F[1] = -m2 * x0_dotdot(t) - Fm(x[0], x[1], 0.0, 0.0)

    # 減衰強制振動モデル
    y = np.dot(M_inv, F) - np.dot((np.dot(C, M_inv)), v) - np.dot((np.dot(K, M_inv)), x)

    return y

def RK4step(v, x, i_f, i_s, t):
    k11 = f(x, v) * dt
    k12 = v * dt
    k21 = f(x + k11 / 2, v + k12 / 2) * dt
    k22 = (v + k12 / 2) * dt
    k31 = f(x + k21 / 2, v + k22 / 2) * dt
    k32 = (v + k22 / 2) * dt
    k41 = f(x + k31, v + k32) * dt
    k42 = (v + k32) * dt
    v += (k11 + 2 * k21 + 2 * k31 + k41) / 6
    x += (k12 + 2 * k22 + 2 * k32 + k42) / 6

    return v, x

def x0(t):
    return Amp * np.sin(2 * np.pi * freq * t)  # 加振器変位

def x0_dot(t):
    return Amp * 2 * np.pi * freq * np.cos(2 * np.pi * freq * t)  # 加振器速度

def x0_dotdot(t):
    # return -Amp * (2 * np.pi * freq)**2 * np.sin(2 * np.pi * freq * t)  # 加振器加速度
    return Amp * np.sin(2 * np.pi * freq * t)  # 加振器加速度

if __name__ == '__main__':
    ''' 固有値解析 '''
    # 質量マトリクスの逆行列を計算
    M_inv = np.linalg.inv(M)

    # 固有値と固有ベクトルを計算
    omega, v = np.linalg.eig(np.dot(M_inv, K))

    # 固有値の順番を昇順にソートして固有振動数[Hz]に変換
    omega_sort = (1 / (2 * np.pi)) * np.sqrt(np.sort(omega))

    # 固有値のソート時のインデックスを取得
    sort_index = np.argsort(omega)

    # 固有値に対応する固有ベクトルをソート
    v_sort = []
    for i in range(len(sort_index)):
        v_sort.append(v.T[sort_index[i]])
    v_sort = np.array(v_sort)

    eigen_value1, eigen_vector1 = omega_sort, v_sort
    print(eigen_value1)

    """数値解析"""
    # 初期化:x[m], v[m/s]
    x = np.zeros(len(M))
    x[0] = +x0_0
    v = np.zeros(len(M))

    t_axis = np.arange(t0, t1, dt)
    x_sol = []
    x_1_sol = []
    v_sol = []

    # 4次のRunge-Kutta法による数値解析(外力成分を離散1D波形から抽出するバージョン）
    iteration = 0
    for t in t_axis:
        # 1自由度目の変位と速度をappend
        x_sol.append(x[0])
        v_sol.append(v[0])
        x_1_sol.append(x[1])
        # イタレーションをモニタ
        if  iteration % 100 == 0:
            print('iteration=', iteration, 'time=', t, x0_dotdot(t))

        # Runge-Kutta
        v, x = RK4step(v, x, i_f, i_s, t)
        iteration += 1

    # 入出力信号のフーリエ変換
    # fft_i = fftpack.fft(x0_dotdot(t))
    # fft_o = fftpack.fft(x_sol)

    # # FRFを計算
    # h_io = (fft_o * fft_i.conjugate()) / (fft_i * fft_i.conjugate())
    # amp = np.sqrt((h_io.real ** 2) + (h_io.imag ** 2))
    # amp = amp / (len(x0_dotdot(t)) / 2)
    # phase = np.arctan2(h_io.imag, h_io.real)
    # phase = np.degrees(phase)
    # freq = np.linspace(0, 1/dt, len(x0_dotdot(t)))

    # プロット
    # フォントの種類とサイズを設定する。
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする。
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右に目盛線を付ける。
    fig = plt.figure(figsize=(10, 6))
    ax2 = fig.add_subplot(111)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')

    # 軸のラベルを設定する。
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Displacement[m]')

    # データプロットの準備とともに、ラベルと線の太さ、凡例の設置を行う。
    x_sol = np.array(x_sol)
    x_1_sol = np.array(x_1_sol)
    # ax1.plot(t_axis, x0_dotdot(t), label='1dof', lw=1, color='red')
    ax2.plot(t_axis, x_sol - x_1_sol, label='z1-z2', lw=1, color='red')
    # ax3.plot(freq, phase, label='1dof', lw=1, color='red')
    # ax4.plot(freq, amp, label='1dof', lw=1, color='red')

    # レイアウト設定
    fig.tight_layout()

    # グラフを表示する。
    plt.show()
    plt.close()
