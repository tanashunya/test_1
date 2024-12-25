import numpy as np
from scipy import fftpack
from scipy.signal import chirp
import matplotlib.pyplot as plt
from RBFN_model import Fm_Linear, Phi_RBFN

# 定数
m1 = 0.1    
m2 = 1.0
c1 = 0.5
c2 = 0.77
c12 = 0.0

# k1 = 200.0
# k2 = 10.0
# k_12 = 10.0

# k1 = 10.0
# k2 = 10.0
# k_12 = 10.0

# k1 = 5.0
# k2 = 0.0
# k_12 = 0.0

# k1 = 0.0
# k2 = 0.0
# k_12 = 0.0


# 室蘭―深い
# k1 = 200.0
# k2 = 10.0
# k_12 = 10.0

# k1 = 10.0
# k2 = 5.0
# k_12 = 5.0

k1 = 20.0
k2 = 10.0
k_12 = 5.0

M = np.array([[m1, 0],
                [0, m2]])
C = np.array([[c1 + c12, -c12],
                [-c12, c2 + c12]])
K = np.array([[k1 + k_12, -k_12],
                [-k_12, k2 + k_12]])

i_f = 0.0
i_s = 0.0
Z = 10000000

# 解析条件
t0 = 0.0
t1 = 1.0
dt = 1e-3
# 初期値
x0_0 = 0.0001

sweep = True
f0 = 3.0
f1 = 0.0001

freq = 1.8
Amp = 0.1

fm_model = Fm_Linear('sum2.csv')
phi_model = Phi_RBFN('Phi_sum_final_model.pth')

def Fm(z1, z2, i_f, i_s):
    # return 0.0
    fm_val = fm_model.interpolate(z1-z2, i_f, i_s)
    # print(z1-z2)
    return fm_val

def Phi(z1, z2, i_f, i_s):
    return phi_model.predict([[z1-z2, i_f, i_s]])

def f(x, v, i_f, i_s, t):
    ''' 運動方程式 '''

    # 質量マトリクスの逆行列
    M_inv = np.linalg.inv(M)

    # 外力ベクトル
    F = np.zeros(len(M))
    F[0] = -m1 * x0_dotdot(t) + Fm(x[0], x[1], i_f, i_s)
    F[1] = -m2 * x0_dotdot(t) - Fm(x[0], x[1], i_f, i_s)

    # 減衰強制振動モデル
    y = np.dot(M_inv, F) - np.dot((np.dot(C, M_inv)), v) - np.dot((np.dot(K, M_inv)), x)

    return y

def RK4step(v, x, i_f, i_s, t):
    k11 = f(x, v, i_f, i_s, t) * dt
    k12 = v * dt
    k21 = f(x + k11 / 2, v + k12 / 2, i_f, i_s, t) * dt
    k22 = (v + k12 / 2) * dt
    k31 = f(x + k21 / 2, v + k22 / 2, i_f, i_s, t) * dt
    k32 = (v + k22 / 2) * dt
    k41 = f(x + k31, v + k32, i_f, i_s, t) * dt
    k42 = (v + k32) * dt
    v += (k11 + 2 * k21 + 2 * k31 + k41) / 6
    x += (k12 + 2 * k22 + 2 * k32 + k42) / 6

    return v, x

def backward_euler_step(phi_current, v, x, i_f, i_s, t):
    # ニュートン法で0にしたい式
    def g(i_s):
        return (Phi(x[0], x[1], i_f, i_s) - phi_current)/dt + Z * i_s

    h           = 1e-10
    NR_conv     = 1e-13
    NR_max_iter = 100

    # ニュートン法開始
    i_s_k = i_s
    for k in range(NR_max_iter):
        G = g(i_s_k)
        dGdi_s = (8*(g(i_s_k + h) - g(i_s_k - h)) - (g(i_s_k + 2*h) - g(i_s_k - 2*h)))/(12*h) #4次中心差分

        i_s_kp1 = i_s_k - G/dGdi_s
        if abs(i_s_kp1 - i_s_k) < NR_conv or k == NR_max_iter - 1:
            # print(f"NR converged at {k+1}th iteration")
            break
        i_s_k = i_s_kp1
    
    i_s_np1 = i_s_kp1
    phi_np1 = Phi(x[0], x[1], i_f, i_s_np1)
    return phi_np1, i_s_np1

def gear_second_step(phi_n, phi_nm1, v, x, i_f, i_s, t):
    
    # ニュートン法で0にしたい式
    def g(i_s):
        return (3*Phi(x[0], x[1], i_f, i_s) - 4*phi_n + phi_nm1)/(2*dt) + Z * i_s

    h           = 1e-10
    NR_conv     = 1e-13
    NR_max_iter = 100

    # ニュートン法開始
    i_s_k = i_s
    for k in range(NR_max_iter):
        G = g(i_s_k)
        # dGdi_s = (g(i_s_k + h) - G)/(h) #1次前進差分
        # dGdi_s = (g(i_s_k + h) - g(i_s_k - h))/(2*h) #2次中心差分
        dGdi_s = (8*(g(i_s_k + h) - g(i_s_k - h)) - (g(i_s_k + 2*h) - g(i_s_k - 2*h)))/(12*h) #4次中心差分

        i_s_kp1 = i_s_k - G/dGdi_s
        if abs(i_s_kp1 - i_s_k) < NR_conv or k == NR_max_iter - 1:
            # print(f"NR converged at {k+1}th iteration")
            break
        i_s_k = i_s_kp1
    
    i_s_np1 = i_s_kp1
    phi_np1 = Phi(x[0], x[1], i_f, i_s_np1)
    return phi_np1, i_s_np1

def x0(t):
    return Amp * np.sin(2 * np.pi * freq * t)  # 加振器変位

def x0_dot(t):
    # return Amp * 2 * np.pi * freq * np.cos(2 * np.pi * freq * t)  # 加振器速度
    return (Amp) * 2 * np.pi * freq * np.cos(2 * np.pi * freq * t)  # 加振器速度

def x0_dotdot(t):
    # return -Amp * (2 * np.pi * freq)**2 * np.sin(2 * np.pi * freq * t)  # 加振器加速度
    if sweep:
        freq_sweep = f0 + (f1 - f0) * t / (t1 - t0)
        return Amp * np.sin(2 * np.pi * freq_sweep * t)  # 加振器加速度
    else:
        return Amp * np.sin(2 * np.pi * freq * t)

def calculate_potential_energy(z1, z2):
    E_mag = 0
    dz = 0.0001
    if z1 >= z2:
        for z in np.arange(0, z1 - z2, dz):
            E_mag -= Fm(z, 0, 0, 0) * dz
    else:
        for z in np.arange(0, z2 - z1, dz):
            E_mag -= Fm(z, 0, 0, 0) * dz

    E = E_mag + 0.5 * k1 * pow(z1, 2) + 0.5 * k2 * pow(z2, 2) + 0.5 * k_12 * pow(z1 - z2, 2)
    return E

if __name__ == '__main__':
    # 数値解析
    # 初期条件
    x = np.zeros(len(M))
    x[0] = +x0_0
    v = np.zeros(len(M))
    phi_current = Phi(x[0], x[1], i_f, i_s)

    t_axis = np.arange(t0, t1, dt)
    x_sol = []
    x_1_sol = []
    v_sol = []

    # 結果を格納する配列
    z1_vals = np.zeros(len(t_axis))
    z2_vals = np.zeros(len(t_axis))
    z1_minus_z2_vals = np.zeros(len(t_axis))
    phi_vals = np.zeros(len(t_axis))
    i_s_vals = np.zeros(len(t_axis))
    Fm_vals = np.zeros(len(t_axis))

    # ループ
    iteration = 0
    for t in t_axis:
        # 1自由度目の変位と速度をappend
        x_sol.append(x[0])
        v_sol.append(v[0])
        x_1_sol.append(x[1])
        # イタレーションをモニタ
        if  iteration % 100 == 0:
            print('iteration=', iteration, 'time=', f'{t:.2f}', 'x=', x[0]-x[1])
        
        # Runge-Kutta
        v, x = RK4step(v, x, i_f, i_s, t)

        if iteration == 0:
            phi_next, i_s = backward_euler_step(phi_current, v, x, i_f, i_s, t)
        else:
            phi_next, i_s = gear_second_step(phi_current, phi_prev, v, x, i_f, i_s, t)
        # 更新
        phi_prev = phi_current
        phi_current = phi_next
        z1_vals[iteration] = x[0]
        z2_vals[iteration] = x[1]
        z1_minus_z2_vals[iteration] = x[0] - x[1]
        phi_vals[iteration] = phi_current
        i_s_vals[iteration] = i_s
        Fm_vals[iteration] = Fm(x[0], x[1], i_f, i_s)

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
    plt.rcParams['font.size'] = 12
    plt.rcParams['font.family'] = 'Times New Roman'
    fig = plt.figure(figsize=(10, 6))

    # 目盛を内側にする。
    # plt.rcParams['xtick.direction'] = 'in'
    # plt.rcParams['ytick.direction'] = 'in'

    # ax2 = fig.add_subplot(111)
    # ax2.yaxis.set_ticks_position('both')
    # ax2.xaxis.set_ticks_position('both')
    # ax2.plot(t_axis, z1_minus_z2_vals, label='z1-z2', lw=1, color='red')
    # ax2.set_xlabel('Time [s]')
    # ax2.set_ylabel('Displacement[m]')

    # データプロットの準備とともに、ラベルと線の太さ、凡例の設置を行う。
    # x_sol = np.array(x_sol)
    # x_1_sol = np.array(x_1_sol)
    # ax1.plot(t_axis, x0_dotdot(t), label='1dof', lw=1, color='red')
    # ax3.plot(freq, phase, label='1dof', lw=1, color='red')
    # ax4.plot(freq, amp, label='1dof', lw=1, color='red')

    # レイアウト設定
    # fig.tight_layout()

    # グラフを表示する。
    # plt.show()
    # plt.close()

    plt.figure()
    plt.plot(t_axis, z1_vals, label='z1')
    plt.plot(t_axis, z2_vals, label='z2')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.legend()
    plt.savefig('png_new/z1_z2_plot.png')

    plt.figure()
    plt.plot(t_axis, z1_minus_z2_vals, label='z1 - z2')
    plt.xlabel('Time [s]')
    plt.ylabel('z1 - z2 [m]')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(10))
    plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.001))
    plt.gca().set_xlim([t_axis[0], t_axis[-1]])  # x軸の範囲を設定
    # plt.gca().set_ylim([-0.02, 0.02])  # y軸の範囲を設定
    plt.savefig('png_new/z1_minus_z2_plot.png')

    plt.figure()
    plt.plot(t_axis, phi_vals, label='Phi')
    plt.xlabel('Time [s]')
    plt.ylabel('Phi [Wb]')
    plt.legend()
    plt.savefig('png_new/phi_plot.png')

    plt.figure()
    plt.plot(t_axis, i_s_vals, label='i_s')
    plt.xlabel('Time [s]')
    plt.ylabel('i_s [A]')
    plt.legend()
    plt.savefig('png_new/i_s_plot.png')

    plt.figure()
    plt.plot(t_axis, Fm_vals, label='Fm')
    plt.xlabel('Time [s]')
    plt.ylabel('Fm [N]')
    plt.legend()
    plt.savefig('png_new/Fm_plot.png')

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

    print("eigen_value1=", eigen_value1)
    print("eigen_vector1=", eigen_vector1)

    # 系のポテンシャルエネルギー構造の算出
    z1_range = np.linspace(-0.02, 0.02, 30)
    z2_range = np.linspace(-0.02, 0.02, 30)
    E = np.zeros((len(z1_range), len(z2_range)))
    for i in range(len(z1_range)):
        for j in range(len(z2_range)):
            E[i, j] = calculate_potential_energy(z1_range[i], z2_range[j])
    # プロット
    plt.figure(figsize=(10, 6))
    contour = plt.contourf(z1_range, z2_range, E, levels=50, cmap='viridis')
    plt.colorbar(contour)
    plt.xlabel('z1 [m]')
    plt.ylabel('z2 [m]')
    # plt.title('Potential Energy Contour')
    plt.savefig('png_new/potential_energy_contour3d.png')


    z1_range = np.linspace(-0.02, 0.02, 100)
    E = np.zeros((len(z1_range)))
    for i in range(len(z1_range)):
        E[i] = calculate_potential_energy(z1_range[i], 0)
    # プロット
    plt.figure(figsize=(10, 6))
    plt.plot(z1_range, E)
    plt.xlabel('z1 [m]')
    plt.ylabel('Potential Energy [J]')
    # plt.title('Potential Energy Contour')
    plt.savefig('png_new/potential_energy_contour2d.png')
