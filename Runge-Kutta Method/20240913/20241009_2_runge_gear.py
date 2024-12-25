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
k1 = 0.0
k2 = 0.0
k_12 = 0.0


# 解析条件
t0 = 0.0
t1 = 10.0
dt = 1e-3

f0 = 4.0
f1 = 4.0
Amp = 0.01

fm_model = Fm_Linear('sum2.csv')

def Fm(z1, z2, i_f, i_s):
    # return 0.0
    fm_val = fm_model.interpolate(z1-z2, i_f, i_s)
    print(z1-z2,fm_val)
    return fm_val

def model_2dof():
    ''' 質量、減衰、剛性マトリクスを定義（2自由度モデル） '''
    
    M = np.array([[m1, 0],
                  [0, m2]])
    C = np.array([[c1 + c12, -c12],
                  [-c12, c2 + c12]])
    K = np.array([[k1 + k_12, -k_12],
                  [-k_12, k2 + k_12]])

    return M, C, K

def f(x, v, f, M, C, K):
    ''' 運動方程式 '''

    # 質量マトリクスの逆行列
    M_inv = np.linalg.inv(M)

    # 外力ベクトル
    F = np.zeros(len(M))
    F[0] = f + Fm(x[0], x[1], 0.0, 0.0)
    F[1] = f - Fm(x[0], x[1], 0.0, 0.0)

    # 減衰強制振動モデル
    y = np.dot(M_inv, F) - np.dot((np.dot(C, M_inv)), v) - np.dot((np.dot(K, M_inv)), x)

    return y


if __name__ == '__main__':
    ''' メイン処理 '''
    # 2自由度振動系の過渡応答解析
    M, C, K = model_2dof()

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



    """数値解析"""
    # 初期化:x0[m], v0[m/s]
    x0 = np.zeros(len(M))
    x0[0] = +0.001
    v0 = np.zeros(len(M))
    x, v = x0, v0

    t_axis = np.arange(t0, t1, dt)
    x_sol = []
    x_1_sol = []
    v_sol = []

    # サインスイープ波形を加振力に設定
    force = Amp * chirp(t_axis, f0=f0, f1=f1, t1=t_axis[-1], method='linear')

    # 4次のRunge-Kutta法による数値解析(外力成分を離散1D波形から抽出するバージョン）
    iteration = 0
    for t in t_axis:
        # 1自由度目の変位と速度をappend
        x_sol.append(x[0])
        v_sol.append(v[0])
        x_1_sol.append(x[1])
        # イタレーションをモニタ
        if  iteration % 100 == 0:
            print('iteration=', iteration, 'time=', t, force[iteration])

        # Runge-Kuttaのメイン
        k11 = f(x, v, force[iteration], M, C, K) * dt
        k12 = v * dt
        k21 = f(x + k11 / 2, v + k12 / 2, force[iteration], M, C, K) * dt
        k22 = (v + k12 / 2) * dt
        k31 = f(x + k21 / 2, v + k22 / 2, force[iteration], M, C, K) * dt
        k32 = (v + k22 / 2) * dt
        k41 = f(x + k31, v + k32, force[iteration], M, C, K) * dt
        k42 = (v + k32) * dt
        v += (k11 + 2 * k21 + 2 * k31 + k41) / 6
        x += (k12 + 2 * k22 + 2 * k32 + k42) / 6
        iteration += 1
    
    force, x_sol, dt, t_axis
    print(eigen_value1)

    # FRFを関数を実行して計算

    # 入出力信号のフーリエ変換
    fft_i = fftpack.fft(force)
    fft_o = fftpack.fft(x_sol)

    # FRFを計算
    h_io = (fft_o * fft_i.conjugate()) / (fft_i * fft_i.conjugate())
    amp = np.sqrt((h_io.real ** 2) + (h_io.imag ** 2))
    amp = amp / (len(force) / 2)
    phase = np.arctan2(h_io.imag, h_io.real)
    phase = np.degrees(phase)
    freq = np.linspace(0, 1/dt, len(force))

    # プロット
    # フォントの種類とサイズを設定する。
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # 目盛を内側にする。
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # グラフの上下左右に目盛線を付ける。
    fig = plt.figure(figsize=(10, 6))
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
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Force[N]')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Displacement[m]')
    ax3.set_xlabel('Excitation frequency [Hz]')
    ax3.set_ylabel('Phase[deg.]')
    ax4.set_xlabel('Excitation frequency [Hz]')
    ax4.set_ylabel('Amplitude[m/N]')

    # スケールの設定をする。
    ax3.set_xlim(0.1, 100)
    ax3.set_yticks(np.arange(-270, 270, 90))
    ax3.set_ylim(-180, 180)
    ax3.set_xscale('log')
    ax4.set_xlim(0.1, 100)
    #ax4.set_ylim(1e-8, 1e-3)
    ax4.loglog()

    # データプロットの準備とともに、ラベルと線の太さ、凡例の設置を行う。
    x_sol = np.array(x_sol)
    x_1_sol = np.array(x_1_sol)
    ax1.plot(t_axis, force, label='1dof', lw=1, color='red')
    ax2.plot(t_axis, x_sol - x_1_sol, label='z1-z2', lw=1, color='red')
    ax3.plot(freq, phase, label='1dof', lw=1, color='red')
    ax4.plot(freq, amp, label='1dof', lw=1, color='red')

    # レイアウト設定
    fig.tight_layout()

    # グラフを表示する。
    plt.show()
    plt.close()
